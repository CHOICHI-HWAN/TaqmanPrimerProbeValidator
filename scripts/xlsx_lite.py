"""Minimal .xlsx reader/writer on the standard library only (zipfile + xml), for the primer / standard-material tables.

Writes one or more sheets of plain text cells (inline strings, bold header row, frozen header, column widths) that Excel,
LibreOffice and Google Sheets open; reads the first sheet of any .xlsx (shared strings, inline strings, numbers, booleans).
No formulas, styles or dates are interpreted -- cell values come back as text, which is what these tables hold.
"""
import io
import re
import zipfile
from xml.etree import ElementTree as ET
from xml.sax.saxutils import escape

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
      "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
      "pr": "http://schemas.openxmlformats.org/package/2006/relationships"}
_ILLEGAL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")


def _col_letter(i):                       # 0 -> A, 26 -> AA
    s = ""
    i += 1
    while i:
        i, r = divmod(i - 1, 26)
        s = chr(65 + r) + s
    return s


def _col_index(ref):                      # "AB12" -> 27
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def _sheet_xml(rows, widths):
    out = ['<?xml version="1.0" encoding="UTF-8" standalone="yes"?>',
           f'<worksheet xmlns="{NS["m"]}"><sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
           'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews><sheetFormatPr defaultRowHeight="15"/>']
    if widths:
        out.append("<cols>" + "".join(f'<col min="{i + 1}" max="{i + 1}" width="{w}" customWidth="1"/>' for i, w in enumerate(widths)) + "</cols>")
    out.append("<sheetData>")
    for ri, row in enumerate(rows, start=1):
        out.append(f'<row r="{ri}">')
        for ci, val in enumerate(row):
            if val is None or val == "":
                continue
            text = _ILLEGAL.sub("", str(val))
            style = ' s="1"' if ri == 1 else ""
            out.append(f'<c r="{_col_letter(ci)}{ri}" t="inlineStr"{style}><is><t xml:space="preserve">{escape(text)}</t></is></c>')
        out.append("</row>")
    out.append("</sheetData></worksheet>")
    return "".join(out)


def write_xlsx(sheets):
    """sheets: [(name, rows, widths)] where rows[0] is the header. -> bytes"""
    names = [re.sub(r"[\[\]\*\?/\\:]", "_", n)[:31] or f"Sheet{i + 1}" for i, (n, _r, _w) in enumerate(sheets)]
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
                   '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
                   '<Default Extension="xml" ContentType="application/xml"/>'
                   '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
                   '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
                   + "".join(f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
                             for i in range(len(sheets))) + "</Types>")
        z.writestr("_rels/.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/></Relationships>')
        z.writestr("xl/workbook.xml", f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><workbook xmlns="{NS["m"]}" xmlns:r="{NS["r"]}"><sheets>'
                   + "".join(f'<sheet name="{escape(n)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>' for i, n in enumerate(names)) + "</sheets></workbook>")
        z.writestr("xl/_rels/workbook.xml.rels", '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
                   '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
                   + "".join(f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
                             for i in range(len(sheets)))
                   + f'<Relationship Id="rId{len(sheets) + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/></Relationships>')
        z.writestr("xl/styles.xml", f'<?xml version="1.0" encoding="UTF-8" standalone="yes"?><styleSheet xmlns="{NS["m"]}">'
                   '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font><font><b/><sz val="11"/><name val="Calibri"/></font></fonts>'
                   '<fills count="2"><fill><patternFill patternType="none"/></fill><fill><patternFill patternType="gray125"/></fill></fills>'
                   '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
                   '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
                   '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
                   '<xf numFmtId="0" fontId="1" fillId="0" borderId="0" xfId="0" applyFont="1"/></cellXfs>'
                   '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles></styleSheet>')
        for i, (_n, rows, widths) in enumerate(sheets):
            z.writestr(f"xl/worksheets/sheet{i + 1}.xml", _sheet_xml(rows, widths))
    return buf.getvalue()


def _num_text(v):
    try:
        f = float(v)
    except ValueError:
        return v
    return str(int(f)) if f == int(f) and "e" not in v.lower() else v


def read_xlsx(data):
    """First sheet of an .xlsx -> (sheet name, rows as lists of strings, padded to the widest row). Raises ValueError."""
    try:
        z = zipfile.ZipFile(io.BytesIO(data))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
    except (zipfile.BadZipFile, KeyError, ET.ParseError):
        raise ValueError("엑셀(.xlsx) 파일이 아니거나 손상되었습니다 (.xls·.csv는 .xlsx로 저장해 주세요)") from None
    first = wb.find("m:sheets/m:sheet", NS)
    if first is None:
        raise ValueError("시트가 없습니다")
    rid = first.get(f'{{{NS["r"]}}}id')
    rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
    target = next((r.get("Target") for r in rels.findall("pr:Relationship", NS) if r.get("Id") == rid), None)
    if target is None:
        raise ValueError("첫 시트를 찾을 수 없습니다")
    path = target.lstrip("/") if target.startswith("/") else "xl/" + target
    shared = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", NS):
            shared.append("".join(t.text or "" for t in si.iter(f'{{{NS["m"]}}}t')))
    sheet = ET.fromstring(z.read(path))
    rows = []
    for row in sheet.findall("m:sheetData/m:row", NS):
        cells = {}
        for c in row.findall("m:c", NS):
            ref, t = c.get("r") or "", c.get("t")
            v = c.find("m:v", NS)
            if t == "inlineStr":
                text = "".join(x.text or "" for x in c.iter(f'{{{NS["m"]}}}t'))
            elif v is None or v.text is None:
                continue
            elif t == "s":
                text = shared[int(v.text)] if int(v.text) < len(shared) else ""
            elif t == "b":
                text = "TRUE" if v.text.strip() == "1" else "FALSE"
            elif t in ("str", "e"):
                text = v.text
            else:
                text = _num_text(v.text)
            if ref:
                cells[_col_index(ref)] = text
        r_idx = int(row.get("r") or len(rows) + 1)
        while len(rows) < r_idx - 1:
            rows.append([])
        rows.append([cells.get(i, "") for i in range(max(cells) + 1)] if cells else [])
    width = max((len(r) for r in rows), default=0)
    return first.get("name") or "Sheet1", [r + [""] * (width - len(r)) for r in rows]
