const state = {
  sets: [], groups: [], activeId: null, activeType: null, genomeSort: { key: "genome_id", dir: 1 },
  problemsOnly: false, msaFilterMode: "all", msaSearch: "",
  groupMode: false, selectedForGroup: new Set(),
  uiState: { notes: {}, set_order: [], group_order: [] },
  lineageSortMode: "name",
};

// ---------------- host scope (사람 / 동물 / 미상 / 전체) ----------------
// Every number and table in the dashboard counts the genomes of ONE host scope. The server decides what the default is: human genomes
// (diagnostic assays are meant for patient samples), or everything when the metadata has no human genome. Laboratory hosts are in no
// scope. The choice is sent as ?host=... on every API read, so each view is simply re-fetched when it changes.
const HOST_SCOPE_LABEL = { human: "사람", animal: "동물", unknown: "미상", all: "전체" };
const HOST_CLASS_LABEL = { human: "사람", animal: "동물", unknown: "미상", lab: "실험 숙주" };
state.hostScope = null;                    // null = the server's primary scope
try { state.hostScope = localStorage.getItem("lasv_host_scope") || null; } catch (e) { /* private mode etc. */ }
(function installHostScopeFetch() {
  const realFetch = window.fetch.bind(window);
  window.fetch = (input, init) => {
    const method = String((init && init.method) || "GET").toUpperCase();
    if (state.hostScope && typeof input === "string" && input.startsWith("/api/") && method === "GET" && !/[?&]host=/.test(input)) {
      input += (input.includes("?") ? "&" : "?") + "host=" + encodeURIComponent(state.hostScope);
    }
    return realFetch(input, init);
  };
})();

function hostScopeSegments() { return Object.entries((state.hostInfo || {}).segments || {}); }
function activeHostScope() {
  const first = hostScopeSegments()[0];
  return state.hostScope || (first && first[1].scope) || "human";
}
function hostCountsText(info, scope) {                       // "S 571 · L 270"
  return Object.entries(info).map(([seg, i]) => `${seg} ${i.scope_counts[scope] ?? 0}`).join(" · ");
}
function renderHostScopeBox() {
  const box = $("#hostScopeBox");
  if (!box) return;
  const segs = state.hostInfo && state.hostInfo.segments;
  if (!segs || !Object.keys(segs).length) { box.hidden = true; return; }
  const info = segs, active = activeHostScope();
  const lab = Object.values(info).reduce((a, i) => a + (i.excluded_lab || 0), 0);
  box.hidden = false;
  box.innerHTML = `
    <div class="host-scope-title">숙주 기준 <span class="note">(커버리지를 세는 게놈)</span></div>
    <div class="host-scope-btns">${["human", "animal", "unknown", "all"].map((sc) =>
      `<button type="button" class="host-scope-btn${sc === active ? " on" : ""}" data-scope="${sc}" title="${HOST_SCOPE_LABEL[sc]}: ${escHtml(hostCountsText(info, sc))}개 게놈">${HOST_SCOPE_LABEL[sc]}</button>`).join("")}</div>
    <div class="host-scope-note">${HOST_SCOPE_LABEL[active]} 기준 ${escHtml(hostCountsText(info, active))}개${lab ? ` · 실험 숙주 ${lab}개 제외` : ""}</div>`;
  box.querySelectorAll(".host-scope-btn").forEach((b) => b.addEventListener("click", () => setHostScope(b.dataset.scope)));
}
async function loadHostScope() {
  try { state.hostInfo = await fetchJSON("/api/host-scope"); } catch (e) { state.hostInfo = null; }
  renderHostScopeBox();
}
async function setHostScope(scope) {
  state.hostScope = scope;
  try { localStorage.setItem("lasv_host_scope", scope); } catch (e) { /* ignore */ }
  await loadHostScope();
  await Promise.all([loadSets(), loadGroups()]);
  const id = state.activeId, type = state.activeType;
  if (type === "set" && id) return selectSet(id);
  if (type === "group" && id) return selectGroup(id);
  if (!$("#overview").hidden) return showOverview();
  // any other page (Venn diagram, re-design, ...) shows numbers of the host scope too: open it again through its own menu button
  const open = MAIN_VIEW_IDS.find((v) => !$(`#${v}`).hidden);
  const btn = open && NAV_VIEW_BTN[open] && document.getElementById(NAV_VIEW_BTN[open]);
  if (btn) btn.click();
  return null;
}
// A one-line statement of what the numbers on the page count (detail / group / overview headers).
function hostScopeBadgeHtml(infoBySeg) {
  const entries = Object.entries(infoBySeg || {});
  if (!entries.length) return "";
  const scope = entries[0][1].scope, other = ["human", "animal", "unknown"].filter((c) => c !== scope && scope !== "all");
  const n = entries.map(([seg, i]) => `${entries.length > 1 ? seg + " " : ""}${i.scope_counts[scope]}`).join(" · ");
  const rest = scope === "all" ? "" : ` (${other.map((c) => `${HOST_CLASS_LABEL[c]} ${entries.reduce((a, [, i]) => a + i.counts[c], 0)}`).join(" · ")}는 이 수치에서 제외)`;
  const lab = entries.reduce((a, [, i]) => a + (i.excluded_lab || 0), 0);
  return `<div class="host-scope-badge" title="왼쪽 사이드바의 '숙주 기준'에서 바꿀 수 있습니다">
    <b>숙주 기준: ${HOST_SCOPE_LABEL[scope]}</b> — 이 화면의 커버리지는 ${scope === "all" ? "사람·동물·미상 전체" : HOST_SCOPE_LABEL[scope] + " 유래"} 게놈 ${escHtml(n)}개 기준입니다${rest}${lab ? `. 실험 숙주 ${lab}개는 어느 범위에도 포함하지 않습니다` : ""}.</div>`;
}

// ---------------------------------------------------------------------------------------
// English / Korean UI toggle. Rather than refactoring every one of the hundreds of
// template-literal strings scattered across this file into abstract i18n keys (a much
// larger, much riskier rewrite), this keeps Korean as the single source of truth in the
// markup and translates in a generic DOM pass after each render: any leaf element (no
// child elements, so replacing its whole text can't clobber nested structure) whose
// rendered text matches a dictionary entry gets swapped, with the original stashed on the
// element so switching back to Korean is a lossless restore rather than a re-render.
// Scope: sidebar navigation, card/section titles, table headers, static form/modal text,
// and tooltips -- the app's structural chrome. Long free-text explanations (the guide
// page's FAQ prose, the per-set automated assessment reasoning) are intentionally left
// Korean-only for now; translating those accurately needs its own dedicated pass, not one
// bundled into a dictionary swap.
const STRINGS_EN = {
  // sidebar navigation
  "전체 파이프라인 목록": "All Pipelines",
  "인터랙티브 검증 대시보드": "Interactive Validation Dashboard",
  "ⓘ in silico 예측 — 실험으로 확인 필요": "ⓘ In silico prediction — confirm experimentally",
  "전체 대시보드": "Overview",
  "+ 새 프라이머 세트 분석": "+ New Primer Set Analysis",
  "📖 파이프라인 개요 / 가이드": "📖 Pipeline Overview / Guide",
  "🧬 리니지별 유전체 다양성": "🧬 Lineage Nucleotide Diversity",
  "🧬 유전자 보존영역": "🧬 Gene Conservation",
  "🗺 메타데이터 / 지도": "🗺 Metadata / Map",
  "⚙ 분석 파라미터 설정": "⚙ Analysis Parameters",
  "📋 변경 이력 (Changelog)": "📋 Changelog",
  "📄 종합 리포트 PDF 생성": "📄 Generate Summary Report (PDF)",
  "데이터 업데이트": "Update Data",
  "분석 다시하기": "Re-run Analysis",
  "📚 문헌 프라이머": "📚 Literature primers",
  "+ 그룹 만들기 (여러 세트 묶기)": "+ Create Group (combine sets)",
  "그룹 이름 (예: VendorB ALL)": "Group name (e.g. VendorB ALL)",
  "그룹으로 묶기": "Create Group",
  "취소": "Cancel",
  "세트 검색 (이름, 벤더)...": "Search sets (name, vendor)...",
  "닫기": "Close",
  "저장": "Save",

  // card / section titles
  "이 도구는 무엇인가요?": "What is this tool?",
  "전체 워크플로우": "Full Workflow",
  "NCBI Primer-BLAST와 무엇이 다른가요?": "How is this different from NCBI Primer-BLAST?",
  "프로세스 호출 다이어그램": "Process Call Diagram",
  "레퍼런스 / 참고문헌 (도구 버전 & 인용)": "References (tool versions & citations)",
  "용어 설명": "Glossary",
  "종합 평가 (참고용 자동 분석 의견)": "Overall Assessment (automated, for reference)",
  "프라이머 / 프로브 서열": "Primer / Probe Sequences",
  "Population 커버리지": "Population Coverage",
  "표준물질(KDCA) 매치": "Standard Material (KDCA) Match",
  "표준물질(KDCA) 상세 결합 정보": "Standard Material (KDCA) Binding Detail",
  "참조유전체 결합위치": "Reference Genome Binding Site",
  "리니지 분포": "Lineage Distribution",
  "게놈별 F / Probe / R 미스매치": "Per-Genome F / Probe / R Mismatch",
  "이차구조 (Self/Cross-dimer · Hairpin)": "Secondary Structure (Self/Cross-dimer · Hairpin)",
  "이차구조 (Self/Cross-dimer · Hairpin, 그룹 전체)": "Secondary Structure (Self/Cross-dimer · Hairpin, whole group)",
  "메모 / 노트": "Notes",
  "메모 / 노트 (그룹 전체)": "Notes (whole group)",
  "결과 해석 가이드": "How to Read These Results",
  "ΔG (결합 자유에너지)란? 어떻게 해석하나요?": "What is ΔG (binding free energy)? How should I read it?",
  "리니지별 유전체 다양성(pixy π) 위에 결합 위치 오버레이": "Binding-site overlay on lineage nucleotide diversity (pixy π)",
  "리니지별 유전체 다양성(pixy π)이란? 어떻게 해석하나요?": "What is lineage nucleotide diversity (pixy π)? How should I read it?",
  "세트별 Population 커버리지": "Population Coverage by Set",
  "합집합 커버리지": "Union Coverage",
  "분석 파라미터": "Analysis Parameters",
  "사용된 프로그램": "Tools Used",
  "이 세트의 검증 기록": "Verification History for This Set",
  "검증 통과 게놈의 리니지 분포": "Lineage Distribution of Validated Genomes",
  "🧪 실험(Wet-lab) 검증": "🧪 Wet-Lab Verification",
  "🧪 실험(Wet-lab) 검증 (그룹 전체)": "🧪 Wet-Lab Verification (whole group)",
  "실험(Wet-lab) 검증 결과 추가 — 그룹 전체": "Add Wet-Lab Verification Result — whole group",
  "전체 게놈 정렬 (F / Probe / R 결합을 한번에 보기)": "Full Genome Alignment (F / Probe / R binding at a glance)",
  "Nextstrain 계통수 대비 커버리지": "Coverage vs. Nextstrain Phylogenetic Tree",
  "이름순": "by name", "건수순": "by count",
  "+ 검증 결과 추가": "+ Add Verification Result",

  // table headers (still-Korean ones; most were already English)
  "평가": "Assessment", "구성 세트": "Member Sets", "국가": "Country", "도구": "Tool",
  "버전": "Version", "상태": "Status", "용도": "Purpose", "인용/출처": "Citation / Source",

  // params page
  "분석 파라미터 설정": "Analysis Parameters",
  "목표 Coverage (%)": "Target Coverage (%)",
  "Primer mismatch 허용 개수 (기본 1)": "Allowed Primer Mismatches (default 1)",
  "Probe mismatch 허용 개수 (기본 0)": "Allowed Probe Mismatches (default 0)",
  "3' 말단 연속 매치 요구 염기수 (기본 2)": "Required 3' Terminal Consecutive Match (bp, default 2)",
  "mfeprimer CPU 코어 수": "mfeprimer CPU cores",

  // new-set form
  "새 프라이머 세트 분석": "New Primer Set Analysis",
  "이름": "Name", "분석 실행": "Run Analysis",
  "config/primers.tsv에 영구 저장 (다음 전체 snakemake 실행에도 포함됨)":
    "Save permanently to config/primers.tsv (included in every future full Snakemake run)",

  // wet-lab verification modal
  "실험(Wet-lab) 검증 결과 추가": "Add Wet-Lab Verification Result",
  "표준물질(양성대조물질)로 실제 프라이머를 테스트한 결과를 인실리코 예측과 비교합니다":
    "Compares an actual wet-lab test against standard material (positive control) with the in-silico prediction",
  "표준물질 선택": "Select Standard Material",
  "검출 여부 (실험 결과)": "Detected? (assay result)",
  "검출 (양성)": "Detected (positive)", "미검출 (음성)": "Not detected (negative)", "판정 보류": "Inconclusive",
  "Ct 값 (선택)": "Ct value (optional)",
  "형광물질(채널)": "Fluorophore (channel)",
  "메모": "Notes",
  "검증 결과 저장": "Save Verification Result",

  // pipeline progress modal
  "파이프라인 실행 중": "Pipeline Running",

  // common tooltips
  "라이트/다크 모드 전환 (Toggle light/dark mode)": "Toggle light/dark mode",
  "이 결과는 서열 기반 컴퓨터 예측입니다. 실제 결합·증폭·진단 성능을 확정하지 않으므로 실험적 검증이 필요합니다.":
    "These results are computational predictions from sequences. They do not establish real binding, amplification or diagnostic performance and must be confirmed experimentally.",
  "Nextstrain에서 최신 서열/메타데이터를 새로 받고 mfeprimer 인덱스를 다시 만듭니다":
    "Fetches the latest sequences/metadata from Nextstrain and rebuilds the mfeprimer index",
  "현재 데이터로 전체 세트 분석을 처음부터 다시 계산합니다": "Recomputes every set from scratch using the current data",
  "모든 세트의 종합 분석 결과를 세트당 1페이지씩 PDF로 정리합니다": "Compiles every set's results into a PDF, one page per set",
  "창을 닫아도 실행은 서버에서 계속됩니다": "Closing this window doesn't stop the run -- it keeps going on the server",
  "이 세트에 저장된 형광물질 값으로 자동 채워집니다 — PCR 플롯의 어느 채널과 비교했는지 기록":
    "Auto-filled from this set's saved fluorophore -- records which channel the PCR plot was compared against",
  "전체화면": "Fullscreen",
  "최소화 / 펼치기": "Minimize / Expand",
  "최대화 (전체화면)": "Maximize (fullscreen)",
  "클릭하면 이 게놈의 결합 상세를 창으로 보여줍니다": "Click to see this genome's binding detail in a window",
  "지정한 참조유전체 위치 구간에 대해, 필터링된 각 게놈이 실제로 갖는 서열을 FASTA로 받습니다":
    "Downloads the actual sequence each filtered genome has over the given reference-genome position range, as FASTA",
  "FASTA로 내보내서 SeaLion(전문 정렬 뷰어)에서 열기": "Export as FASTA to open in SeaLion (a dedicated alignment viewer)",
  "검색어, 필터, 정렬을 모두 초기 상태로": "Reset search, filters, and sorting",
  "정렬 방식 전환": "Toggle sort order",
  "이름 변경": "Rename", "그룹 이름 변경": "Rename group",
  "노트 있음": "Has a note",
  "이미지 삭제": "Delete image", "검증 기록 삭제": "Delete this record",
};

let __i18nLang = null;
function currentLang() {
  if (__i18nLang) return __i18nLang;
  try { __i18nLang = localStorage.getItem("lasv_lang") || "en"; } catch { __i18nLang = "en"; }
  return __i18nLang;
}
function translateLeaf(elm, attr) {
  const key = attr ? elm.getAttribute(attr) : elm.textContent;
  if (key == null) return;
  const dataKey = `i18n${attr ? "Attr" + attr[0].toUpperCase() + attr.slice(1) : "Text"}`;
  if (currentLang() === "en") {
    const en = STRINGS_EN[elm.dataset[dataKey] || key];
    if (en) {
      if (!elm.dataset[dataKey]) elm.dataset[dataKey] = key;
      if (attr) elm.setAttribute(attr, en); else elm.textContent = en;
    }
  } else if (elm.dataset[dataKey]) {
    if (attr) elm.setAttribute(attr, elm.dataset[dataKey]); else elm.textContent = elm.dataset[dataKey];
  }
}
// Deliberately NOT a blanket "div, span" scan of the whole subtree -- pages like the
// per-genome mismatch table or the full-genome alignment view can hold thousands of data
// cells, and almost none of them are chrome text this dictionary would ever match, so
// scanning them on every render is pure waste. This selector list is every element type
// STRINGS_EN actually has entries for.
const TRANSLATE_TEXT_SELECTOR =
  "h1, h2, h3, button, label, option, summary, th, .note, .modal-title, .modal-sub, .section-title, #lineageSortLabel";

// A few h3 headers mix a title with a nested control (a sort-toggle or "+ add" button) --
// translateLeaf()'s children.length===0 guard correctly skips those (never risk mangling
// the nested button), so this handles just that one shape instead: find the leading text
// node and translate only it, leaving every child element untouched.
function translateLeadingText(elm) {
  const node = [...elm.childNodes].find((n) => n.nodeType === 3 && n.textContent.trim());
  if (!node) return;
  const orig = elm.dataset.i18nLead || node.textContent.trim();
  if (currentLang() === "en") {
    const en = STRINGS_EN[orig];
    if (en) {
      if (!elm.dataset.i18nLead) elm.dataset.i18nLead = orig;
      node.textContent = node.textContent.replace(orig, en);
    }
  } else if (elm.dataset.i18nLead) {
    node.textContent = node.textContent.replace(node.textContent.trim(), elm.dataset.i18nLead);
  }
}

function translateNode(root) {
  root.querySelectorAll(TRANSLATE_TEXT_SELECTOR).forEach((elm) => {
    if (!elm.textContent.trim()) return;
    if (elm.children.length === 0) translateLeaf(elm);
    else if ((elm.tagName === "H3" || elm.tagName === "LABEL") && elm.children.length === 1) translateLeadingText(elm);
  });
  root.querySelectorAll("[title]").forEach((elm) => translateLeaf(elm, "title"));
  root.querySelectorAll("[placeholder]").forEach((elm) => translateLeaf(elm, "placeholder"));
}
function setLang(lang) {
  __i18nLang = lang;
  try { localStorage.setItem("lasv_lang", lang); } catch { /* private mode etc. */ }
  document.documentElement.lang = lang;
  translateNode(document.body);
  const btn = document.querySelector("#langToggleBtn");
  if (btn) btn.textContent = lang === "en" ? "KO" : "EN";
}

function applyCustomOrder(items, order, idFn) {
  // Sort `items` (already in the default order, e.g. alphabetical) by their id's position
  // in `order` (a user-dragged id list) -- items missing from `order` keep their relative
  // default order, appended after every item that IS in `order`.
  if (!order || !order.length) return items;
  const rank = new Map(order.map((id, i) => [id, i]));
  const known = items.filter((it) => rank.has(idFn(it)));
  const unknown = items.filter((it) => !rank.has(idFn(it)));
  known.sort((a, b) => rank.get(idFn(a)) - rank.get(idFn(b)));
  return known.concat(unknown);
}

const $ = (sel) => document.querySelector(sel);
const el = (tag, cls, html) => {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html !== undefined) e.innerHTML = html;
  return e;
};

async function fetchJSON(url, opts) {
  const res = await fetch(url, opts);
  if (!res.ok) throw new Error(`${url} -> ${res.status}`);
  return res.json();
}

function fmtPct(v) { return `${Number(v).toFixed(1)}%`; }

// ---------------- modal dialogs (replace the browser's native alert/confirm/prompt) ----------------
// One small promise-based dialog on top of the same .modal-overlay/.modal-box look the pipeline
// and wet-lab modals use: Esc / backdrop click = cancel, Enter = the primary button.
function uiDialog({ title, message, html, input, defaultValue, placeholder, buttons }) {
  return new Promise((resolve) => {
    const prevFocus = document.activeElement;
    const overlay = el("div", "modal-overlay ui-dialog-overlay");
    const btns = buttons.map((b, i) =>
      `<button type="button" class="btn ${b.primary ? (b.danger ? "btn-danger" : "btn-primary") : "btn-secondary"}" data-i="${i}">${escHtml(b.label)}</button>`).join("");
    overlay.innerHTML = `<div class="modal-box ui-dialog" role="dialog" aria-modal="true">
      <div class="modal-head"><div class="modal-title">${escHtml(title || "")}</div></div>
      <div class="ui-dialog-body">
        ${html ? `<div class="ui-dialog-msg ui-dialog-html">${html}</div>` : `<div class="ui-dialog-msg">${escHtml(message || "")}</div>`}
        ${input ? `<input type="text" class="ui-dialog-input" value="${escHtml(defaultValue || "")}" placeholder="${escHtml(placeholder || "")}" spellcheck="false" autocomplete="off">` : ""}
      </div>
      <div class="ui-dialog-actions">${btns}</div></div>`;
    modalHost().appendChild(overlay);
    const inputEl = overlay.querySelector(".ui-dialog-input");
    let done = false;
    const finish = (b) => {
      if (done) return;
      done = true;
      document.removeEventListener("keydown", onKey, true);
      overlay.remove();
      if (prevFocus && prevFocus.focus) { try { prevFocus.focus(); } catch { /* element gone */ } }
      resolve(b.value(inputEl ? inputEl.value : undefined));
    };
    const cancelBtn = buttons.find((b) => b.cancel) || buttons[buttons.length - 1];
    const primaryBtn = buttons.find((b) => b.primary) || buttons[0];
    const onKey = (ev) => {
      if (ev.key === "Escape") { ev.preventDefault(); ev.stopPropagation(); finish(cancelBtn); }
      else if (ev.key === "Enter" && !(document.activeElement && document.activeElement.tagName === "BUTTON")) {
        ev.preventDefault(); ev.stopPropagation(); finish(primaryBtn);
      }
    };
    document.addEventListener("keydown", onKey, true);
    overlay.addEventListener("mousedown", (ev) => { if (ev.target === overlay) finish(cancelBtn); });
    overlay.querySelectorAll("button[data-i]").forEach((b) => b.addEventListener("click", () => finish(buttons[+b.dataset.i])));
    (inputEl || overlay.querySelector(".btn-primary, .btn-danger") || overlay.querySelector("button")).focus();
    if (inputEl) inputEl.select();
  });
}

// Large read-only window for a drill-down (a genome's binding detail, ...) -- Esc, the X button or a
// click on the backdrop closes it. Replaces the old "expand a row inline below the table" pattern.
// Everything that is not in the page's own colour/text styles is dropped: only the properties that decide how an element LOOKS are copied
// into the clone that is rasterised (keeps the SVG small; layout is recomputed by the browser from the same widths).
const DM_STYLE_PROPS = ["color", "background-color", "background-image", "border-top", "border-right", "border-bottom", "border-left", "border-radius",
  "box-shadow", "font-family", "font-size", "font-weight", "font-style", "line-height", "letter-spacing", "text-align", "text-decoration",
  "text-transform", "white-space", "word-break", "overflow-wrap", "vertical-align", "padding-top", "padding-right", "padding-bottom", "padding-left",
  "margin-top", "margin-right", "margin-bottom", "margin-left", "display", "flex-direction", "flex-wrap", "flex-grow", "flex-shrink", "flex-basis",
  "align-items", "justify-content", "gap", "grid-template-columns", "grid-template-rows", "position", "top", "left", "right", "bottom", "opacity",
  "box-sizing", "table-layout", "border-collapse", "border-spacing", "min-width", "max-width", "list-style-type", "text-overflow", "zoom", "fill", "stroke", "stroke-width"];

async function domToPngBlob(root, maxScale = 2) {
  const all = [root, ...root.querySelectorAll("*")];
  if (all.length > 9000) throw new Error("내용이 너무 커서(요소 9000개 초과) 이미지로 만들 수 없습니다. 텍스트 복사를 쓰세요.");
  const w = Math.ceil(root.offsetWidth), h = Math.ceil(root.scrollHeight);
  const clone = root.cloneNode(true);
  const dst = [clone, ...clone.querySelectorAll("*")];
  all.forEach((src, i) => {
    const d = dst[i], cs = getComputedStyle(src);
    let css = "";
    for (const prop of DM_STYLE_PROPS) css += `${prop}:${cs.getPropertyValue(prop)};`;
    // a scrolling area is drawn at its full size, nothing clipped
    if (/(auto|scroll|hidden)/.test(cs.overflowY + cs.overflowX) && (src.scrollHeight > src.clientHeight + 1 || src.scrollWidth > src.clientWidth + 1)) {
      css += `overflow:visible;max-height:none;height:${src.scrollHeight}px;width:${Math.max(src.scrollWidth, src.clientWidth)}px;`;
    }
    d.setAttribute("style", css);
    if (src.tagName === "INPUT") d.setAttribute("value", src.value);
  });
  clone.querySelectorAll(".dm-tools-only, .detail-modal-tools").forEach((x) => x.remove());
  clone.style.cssText += `;width:${w}px;max-height:none;height:auto;overflow:visible;`;
  const bg = getComputedStyle(document.body).backgroundColor || "#fff";
  const xml = new XMLSerializer().serializeToString(clone);
  const svg = `<svg xmlns="http://www.w3.org/2000/svg" width="${w}" height="${h}"><foreignObject x="0" y="0" width="100%" height="100%">${xml}</foreignObject></svg>`;
  const img = new Image();
  await new Promise((res, rej) => { img.onload = res; img.onerror = () => rej(new Error("이미지를 만들지 못했습니다.")); img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(svg); });
  const scale = Math.min(maxScale, 16000 / w, 16000 / h, Math.sqrt(80e6 / (w * h)));
  const canvas = document.createElement("canvas");
  canvas.width = Math.max(1, Math.floor(w * scale)); canvas.height = Math.max(1, Math.floor(h * scale));
  const ctx = canvas.getContext("2d");
  ctx.fillStyle = bg; ctx.fillRect(0, 0, canvas.width, canvas.height);
  ctx.scale(scale, scale); ctx.drawImage(img, 0, 0);
  return new Promise((res, rej) => canvas.toBlob((blob) => (blob ? res(blob) : rej(new Error("이미지가 너무 큽니다."))), "image/png"));
}

function downloadBlob(blob, name) {
  const url = URL.createObjectURL(blob), a = document.createElement("a");
  a.href = url; a.download = name; document.body.appendChild(a); a.click(); a.remove();
  setTimeout(() => URL.revokeObjectURL(url), 2000);
}

// compact: only minimise / maximise / close (for views that bring their own export buttons, e.g. the MSA modal)
// 세트 상세 화면을 읽는 순서대로 다시 배열하고 위쪽에 구역 이동 막대를 붙입니다: 요약 → 결합 근거(게놈별 표 + MSA) → 표준물질 → 다양성·계통 → 실험·메모 → Probe 설계.
// 카드는 옮기기만 하므로 이미 연결된 기능(MSA↔표 연동 등)은 그대로입니다.
function organizeSetDetail(detail) {
  if (detail.querySelector(".detail-tabs") || !detail.querySelector(".detail-header")) return;
  const kids = [...detail.children];
  const heading = (k) => (k.querySelector(":scope > h3, :scope > .msa-card-header h3, :scope h3") || {}).textContent || "";
  const has = (k, re) => re.test(heading(k));
  const rules = [
    ["summary", "요약", (k) => k.classList.contains("assessment-card") || k.classList.contains("card-grid")],
    ["evidence", "결합 근거", (k) => has(k, /게놈별 F \/ Probe \/ R 미스매치/) || k.classList.contains("msa-card")],
    ["standard", "표준물질", (k) => has(k, /표준물질\(KDCA\) 상세 결합/)],
    ["context", "다양성·계통", (k) => has(k, /pixy/) || k.classList.contains("tree-card")],
    ["wet", "실험·메모", (k) => has(k, /실험\(Wet-lab\)/) || has(k, /메모 \/ 노트/)],
    ["design", "Probe 설계", (k) => k.classList.contains("probe-design-card")],
  ];
  const buckets = rules.map(() => []);
  for (const k of kids) {
    for (let i = 0; i < rules.length; i++) if (rules[i][2](k)) { buckets[i].push(k); break; }
  }
  const evidence = buckets[1];
  evidence.sort((a, b) => (a.classList.contains("msa-card") ? 1 : 0) - (b.classList.contains("msa-card") ? 1 : 0));
  const first = kids.find((k) => buckets.some((b) => b.includes(k)));
  if (!first) return;
  const wrap = el("div", "detail-sections");
  const nav = el("div", "detail-tabs");
  nav.innerHTML = rules.map((r, i) => (buckets[i].length ? `<button type="button" class="detail-tab" data-sect="${r[0]}">${r[1]}</button>` : "")).join("");
  wrap.appendChild(nav);
  rules.forEach((r, i) => {
    if (!buckets[i].length) return;
    const isDesign = r[0] === "design";
    const sect = el(isDesign ? "details" : "section", "detail-sect");
    sect.id = "sect-" + r[0];
    if (isDesign) {
      sect.innerHTML = `<summary class="sect-title">${r[1]} <span class="note">(새 Probe 후보 설계·검증 — 위 "Probe 설계" 탭으로 바로 이동)</span></summary>`;
      try { sect.open = localStorage.getItem("lasv_sect_design_open") !== "0"; } catch (e) { /* ignore */ }
      sect.addEventListener("toggle", () => { try { localStorage.setItem("lasv_sect_design_open", sect.open ? "1" : "0"); } catch (e) { /* ignore */ } });
    } else {
      sect.innerHTML = `<h2 class="sect-title">${r[1]}</h2>`;
    }
    buckets[i].forEach((k) => sect.appendChild(k));
    wrap.appendChild(sect);
  });
  const anchor = kids.find((k) => k.parentNode === detail && !buckets.some((bk) => bk.includes(k)) && kids.indexOf(k) > kids.indexOf(first));
  detail.insertBefore(wrap, anchor || null);
  nav.addEventListener("click", (e) => {
    const b = e.target.closest(".detail-tab");
    if (!b) return;
    const t = detail.querySelector("#sect-" + b.dataset.sect);
    if (t && t.tagName === "DETAILS") t.open = true;
    if (t) t.scrollIntoView({ behavior: "smooth", block: "start" });
  });
  const main = document.querySelector(".main");
  const spy = () => {
    if (!document.body.contains(nav)) { main.removeEventListener("scroll", spy); return; }
    const top = nav.getBoundingClientRect().bottom + 8;
    let cur = null;
    wrap.querySelectorAll(".detail-sect").forEach((x) => { if (x.getBoundingClientRect().top <= top + 40) cur = x.id.replace("sect-", ""); });
    nav.querySelectorAll(".detail-tab").forEach((b) => b.classList.toggle("on", b.dataset.sect === cur));
  };
  main.addEventListener("scroll", spy, { passive: true });
  spy();
}

(function wireFinalResultBoxes() {                         // whenever a detail page is rendered, fill its final-result box
  new MutationObserver(() => {
    const box = document.querySelector("#detail .final-result-box:empty");
    if (box) loadFinalResults().then(() => initFinalResult(box.parentElement));
    const d = document.getElementById("detail");
    if (d.querySelector(".detail-header") && d.querySelector(".msa-card") && !d.querySelector(".detail-tabs")) organizeSetDetail(d);
  }).observe(document.getElementById("detail"), { childList: true });
})();

function openDetailModal({ title, subtitle, html, compact }) {
  const prevFocus = document.activeElement;
  const overlay = el("div", "modal-overlay detail-modal-overlay");
  overlay.innerHTML = `<div class="modal-box detail-modal" role="dialog" aria-modal="true">
    <div class="modal-head"><div style="min-width:0;"><div class="modal-title">${title}</div>${subtitle ? `<div class="modal-sub">${subtitle}</div>` : ""}</div>
      <div class="detail-modal-tools">
        ${compact ? "" : `<button type="button" class="btn btn-secondary dm-copy dm-hide-min" title="창의 내용을 텍스트로 복사합니다">&#128203; 복사</button>
        <button type="button" class="btn btn-secondary dm-img dm-hide-min" title="창 전체(스크롤 밖 부분 포함)를 PNG 이미지로 저장합니다">&#128247; 화면 저장</button>
        <button type="button" class="btn btn-secondary dm-copyimg dm-hide-min" title="창 전체를 이미지로 클립보드에 복사합니다 (HTTPS/localhost에서만, 아니면 저장)">&#128444; 이미지 복사</button>
        <button type="button" class="btn btn-secondary dm-fit dm-hide-min on" title="내용이 창보다 넓으면 창을 넓히고, 그래도 넘치면 축소해서 맞춥니다 (끄면 가로 스크롤)">&#8660; 맞춤</button>`}
        <button type="button" class="btn btn-secondary dm-min" title="최소화 (왼쪽 아래로 접기)">&#9473;</button>
        <button type="button" class="btn btn-secondary dm-max dm-hide-min" title="최대화 / 원래 크기">&#9974;</button>
        <button type="button" class="btn btn-secondary detail-modal-close" title="닫기 (Esc)">닫기</button></div></div>
    <div class="detail-modal-body">${html}</div></div>`;
  modalHost().appendChild(overlay);
  const box = overlay.querySelector(".detail-modal"), body = overlay.querySelector(".detail-modal-body");
  let fitOn = !compact;
  const fit = () => {                                      // the window follows its content: wider if needed, then scaled down
    if (overlay.classList.contains("minimized")) return;
    body.style.zoom = ""; box.style.width = "";
    if (!fitOn) return;
    const maxW = window.innerWidth - 56;
    if (!overlay.classList.contains("maximized")) {
      const need = body.scrollWidth + 40;
      if (need > box.getBoundingClientRect().width) box.style.width = `${Math.min(need, maxW)}px`;
    }
    if (body.scrollWidth > body.clientWidth + 2) body.style.zoom = String(Math.max(0.5, body.clientWidth / body.scrollWidth).toFixed(3));
  };
  const close = () => {
    document.removeEventListener("keydown", onKey, true);
    window.removeEventListener("resize", fit);
    overlay.remove();
    if (prevFocus && prevFocus.focus) { try { prevFocus.focus(); } catch { /* element gone */ } }
  };
  const onKey = (ev) => { if (ev.key === "Escape") { ev.preventDefault(); ev.stopPropagation(); close(); } };
  document.addEventListener("keydown", onKey, true);
  window.addEventListener("resize", fit);
  overlay.addEventListener("mousedown", (ev) => { if (ev.target === overlay) close(); });
  overlay.querySelector(".detail-modal-close").addEventListener("click", close);
  overlay.querySelector(".detail-modal-close").focus();
  const flash = (btn, text) => { const o = btn.innerHTML; btn.textContent = text; setTimeout(() => { btn.innerHTML = o; }, 1400); };
  const tool = (cls, fn) => { const b = overlay.querySelector(cls); if (b) b.addEventListener("click", (e) => { e.stopPropagation(); fn(b); }); };
  tool(".dm-copy", (b) => copyTextToClipboard(body.innerText).then(() => flash(b, "✓ 복사됨")).catch(() => uiAlert("복사하지 못했습니다.")));
  const fname = () => `${(overlay.querySelector(".modal-title").textContent || "detail").trim().replace(/[^\w가-힣.-]+/g, "_").slice(0, 60)}.png`;
  tool(".dm-img", async (b) => { flash(b, "저장 중..."); try { downloadBlob(await domToPngBlob(box), fname()); } catch (err) { uiAlert(err.message); } });
  tool(".dm-copyimg", async (b) => {
    try {
      const blob = await domToPngBlob(box);
      if (navigator.clipboard && window.ClipboardItem && window.isSecureContext) { await navigator.clipboard.write([new ClipboardItem({ "image/png": blob })]); flash(b, "✓ 복사됨"); }
      else { downloadBlob(blob, fname()); uiAlert("이 주소(http)에서는 브라우저가 이미지 복사를 허용하지 않아 PNG 파일로 저장했습니다. (localhost 또는 HTTPS에서는 바로 복사됩니다.)"); }
    } catch (err) { uiAlert(err.message); }
  });
  tool(".dm-fit", (b) => { fitOn = !fitOn; b.classList.toggle("on", fitOn); fit(); });
  tool(".dm-max", () => { overlay.classList.toggle("maximized"); fit(); });
  const minBtn = overlay.querySelector(".dm-min");
  const setMin = (on) => {
    overlay.classList.toggle("minimized", on);
    minBtn.innerHTML = on ? "&#9633;" : "&#9473;"; minBtn.title = on ? "복원" : "최소화 (왼쪽 아래로 접기)";
    if (on) {                                              // dock next to the other minimised windows
      const idx = [...document.querySelectorAll(".detail-modal-overlay.minimized")].indexOf(overlay);
      overlay.style.left = `${12 + Math.max(0, idx) * 352}px`;
    } else overlay.style.left = "";
    fit();
  };
  minBtn.addEventListener("click", (e) => { e.stopPropagation(); setMin(!overlay.classList.contains("minimized")); });
  box.addEventListener("click", (e) => { if (overlay.classList.contains("minimized") && !e.target.closest("button")) setMin(false); });
  requestAnimationFrame(fit);
  return { close, overlay, fit };
}

// ---------------- toasts (작업 진행/완료 알림) ----------------
// toast(msg, {type: "info"|"run"|"ok"|"err", id, ttl}) -- same id updates the existing toast in place; "run" is sticky with a spinner.
function toast(msg, opts = {}) {
  let box = document.getElementById("toastBox");
  if (!box) { box = document.createElement("div"); box.id = "toastBox"; box.setAttribute("aria-live", "polite"); document.body.appendChild(box); }
  const type = opts.type || "info", id = opts.id || "t" + Math.random().toString(36).slice(2);
  let el = [...box.children].find((x) => x.dataset.id === id);
  if (!el) {
    el = document.createElement("div"); el.dataset.id = id;
    el.innerHTML = '<span class="toast-ic"></span><span class="toast-msg"></span><button class="toast-x" type="button" aria-label="닫기">×</button>';
    el.querySelector(".toast-x").addEventListener("click", () => toastClear(id));
    box.appendChild(el);
  }
  el.className = `toast toast-${type}`;
  el.querySelector(".toast-ic").textContent = type === "ok" ? "✓" : type === "err" ? "!" : type === "run" ? "" : "i";
  el.querySelector(".toast-msg").textContent = msg;
  clearTimeout(el._t);
  const ttl = opts.ttl ?? (type === "run" ? 0 : type === "err" ? 12000 : 4500);
  if (ttl) el._t = setTimeout(() => toastClear(id), ttl);
  return id;
}
function toastClear(id) {
  const el = [...(document.getElementById("toastBox")?.children || [])].find((x) => x.dataset.id === id);
  if (el) { clearTimeout(el._t); el.classList.add("out"); setTimeout(() => el.remove(), 200); }
}
// Shows a running toast for a job poll loop; call progress(job) each tick and done(ok, msg) at the end.
function jobToast(label) {
  const id = "job-" + label + "-" + Date.now();
  toast(`${label} 시작...`, { type: "run", id });
  return {
    progress: (j) => toast(`${label}: ${j.progress || "실행 중"}${j.elapsed_sec != null ? ` (${fmtElapsedShort(j.elapsed_sec)})` : ""}`, { type: "run", id }),
    done: (ok, msg) => toast(msg || (ok ? `${label} 완료` : `${label} 실패`), { type: ok ? "ok" : "err", id }),
  };
}
function fmtElapsedShort(s) { s = Math.round(Number(s) || 0); return s >= 60 ? `${Math.floor(s / 60)}분 ${s % 60}초` : `${s}초`; }
// Jobs started elsewhere (another computer) show up as one toast too.
(function watchActiveJobs() {
  const KIND = { redesign: "재설계", report: "리포트 생성", pipeline: "파이프라인", analyze: "분석" };
  let busyLocal = 0;
  const tick = async () => {
    try {
      const jobs = await fetchJSON("/api/jobs/active");
      if (Array.isArray(jobs) && jobs.length) toast("서버에서 실행 중: " + jobs.map((j) => `${KIND[j.kind] || j.kind || "작업"}${j.progress ? " – " + j.progress : ""}`).join(" · "), { type: "run", id: "active-jobs" });
      else toastClear("active-jobs");
    } catch (e) { /* offline */ }
  };
  setTimeout(() => { tick(); setInterval(tick, 8000); }, 3000);
})();

function uiAlert(message, opts = {}) {
  return uiDialog({ title: opts.title || "알림", message,
    buttons: [{ label: "확인", primary: true, cancel: true, value: () => undefined }] });
}

function uiConfirm(message, opts = {}) {
  const destructive = opts.danger ?? /삭제|사라집니다/.test(message);
  return uiDialog({ title: opts.title || (destructive ? "삭제 확인" : "확인"), message, html: opts.html,
    buttons: [
      { label: opts.cancelLabel || "취소", cancel: true, value: () => false },
      { label: opts.okLabel || (destructive ? "삭제" : "확인"), primary: true, danger: destructive, value: () => true },
    ] });
}

function uiPrompt(message, defaultValue, opts = {}) {
  return uiDialog({ title: opts.title || "이름 입력", message, input: true, defaultValue,
    buttons: [
      { label: "취소", cancel: true, value: () => null },
      { label: "확인", primary: true, value: (v) => v },
    ] });
}

function escHtml(s) {
  return String(s ?? "").replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));
}

// ---------------- sidebar ----------------

// ---- sidebar tree: group > its sets > (variants as a badge that unfolds) ---------------------------------------------------------
// A set is shown under ONE group -- the smallest group that contains it, so an "ALL" group does not repeat every set -- and sets that
// belong to no group are listed by vendor. A variant ("<parent name>-<n>", see set_parents) is folded into its parent's row.
state.treeClosed = new Set();            // keys "g:<group id>" of groups the user folded
state.treeClosedVariants = new Set();    // keys "v:<set id>" of parents whose variants the user folded (they are shown under their parent by default)
try {
  state.treeClosed = new Set(JSON.parse(localStorage.getItem("lasv_tree_closed") || "[]"));
  state.treeClosedVariants = new Set(JSON.parse(localStorage.getItem("lasv_tree_closed_variants") || "[]"));
} catch (e) { /* private mode */ }
function saveTreeState() {
  try {
    localStorage.setItem("lasv_tree_closed", JSON.stringify([...state.treeClosed]));
    localStorage.setItem("lasv_tree_closed_variants", JSON.stringify([...state.treeClosedVariants]));
  } catch (e) { /* ignore */ }
}
function vendorLabel(source) {
  if (!source) return "출처 없음";
  return /^문헌/.test(source) ? "문헌" : source;
}
function primaryGroupOf(groups) {        // set id -> the smallest group containing it (ties: lowest id)
  const map = new Map();
  [...groups].sort((x, y) => (x.members || []).length - (y.members || []).length || x.group_id.localeCompare(y.group_id)).forEach((g) => {
    for (const id of g.members || []) if (!map.has(id)) map.set(id, g.group_id);
  });
  return map;
}
// the order the overview uses too: [{type:"group", g, sets:[{s, variants}]}, {type:"vendor", name, sets:[...]}] per segment
function buildSetTree(seg, q) {
  const needle = (q || "").toLowerCase();
  const match = (s) => !needle || (s.display_name || s.set_id).toLowerCase().includes(needle) || s.set_id.toLowerCase().includes(needle) || (s.source || "").toLowerCase().includes(needle);
  const all = state.sets.filter((s) => s.segment === seg);
  const ids = new Set(all.map((s) => s.set_id));
  const variantsOf = new Map();
  for (const s of all) if (s.parent_id && ids.has(s.parent_id)) (variantsOf.get(s.parent_id) || variantsOf.set(s.parent_id, []).get(s.parent_id)).push(s);
  const byName = (x, y) => (x.display_name || x.set_id).localeCompare(y.display_name || y.set_id, undefined, { numeric: true });
  const order = state.uiState.set_order;
  const tops = all.filter((s) => !(s.parent_id && ids.has(s.parent_id)));
  const entry = (s) => {
    const vars = (variantsOf.get(s.set_id) || []).slice().sort(byName);
    const shownVars = vars.filter(match);
    if (!match(s) && !shownVars.length) return null;
    return { s, variants: vars, shownVariants: needle ? shownVars : vars };
  };
  const groups = state.groups.filter((g) => g.segment === seg);
  const prim = primaryGroupOf(groups);
  const sortSets = (arr) => applyCustomOrder(arr.slice().sort(byName), order, (s) => s.set_id);
  const out = [];
  const orderedGroups = applyCustomOrder(groups.slice().sort((x, y) => x.label.localeCompare(y.label)), state.uiState.group_order, (g) => g.group_id);
  for (const g of orderedGroups) {
    const kids = sortSets(tops.filter((s) => prim.get(s.set_id) === g.group_id)).map(entry).filter(Boolean);
    const others = groups.filter((x) => x.group_id !== g.group_id && (g.members || []).some((m) => prim.get(m) === x.group_id && x.group_id !== g.group_id));
    const childGroups = [...new Set((g.members || []).map((m) => prim.get(m)).filter((gid) => gid && gid !== g.group_id))].map((gid) => (groups.find((x) => x.group_id === gid) || {}).label).filter(Boolean);
    if (needle && !g.label.toLowerCase().includes(needle) && !kids.length) continue;
    out.push({ type: "group", g, kids, childGroups });
  }
  const free = tops.filter((s) => !prim.has(s.set_id));
  const buckets = {};
  for (const s of free) (buckets[vendorLabel(s.source)] = buckets[vendorLabel(s.source)] || []).push(s);
  for (const name of Object.keys(buckets).sort((x, y) => x.localeCompare(y))) {
    const sets = sortSets(buckets[name]).map(entry).filter(Boolean);
    if (sets.length) out.push({ type: "vendor", name, sets });
  }
  return out;
}

function renderSidebar() {
  const q = ($("#searchBox").value || "").toLowerCase();
  const list = $("#setList");
  list.innerHTML = "";
  const dragging = { id: null, kind: null, container: null };
  const searching = !!q;

  // persist a drag: the dropped item takes the slot of the one it replaced, every other id keeps its place in the saved order
  const commitOrder = (kind, containerIds, newIds) => {
    const key = kind === "set" ? "set_order" : "group_order";
    const base = (state.uiState[key] || []).slice();
    for (const id of containerIds) if (!base.includes(id)) base.push(id);
    const inContainer = new Set(containerIds);
    const slots = base.map((id, i) => (inContainer.has(id) ? i : -1)).filter((i) => i >= 0);
    newIds.forEach((id, k) => { base[slots[k]] = id; });
    state.uiState[key] = base;
    saveOrder(kind, base);
    renderSidebar();
  };
  const makeDraggable = (el2, id, kind, containerKey, idsInOrder) => {
    el2.draggable = true;
    el2.addEventListener("dragstart", (e) => {
      e.stopPropagation();
      dragging.id = id; dragging.kind = kind; dragging.container = containerKey;
      e.dataTransfer.effectAllowed = "move";
    });
    el2.addEventListener("dragover", (e) => {
      if (dragging.kind !== kind || dragging.container !== containerKey || dragging.id === id) return;
      e.preventDefault();
      el2.classList.add("drag-over");
    });
    el2.addEventListener("dragleave", () => el2.classList.remove("drag-over"));
    el2.addEventListener("drop", (e) => {
      e.preventDefault(); e.stopPropagation();
      el2.classList.remove("drag-over");
      if (dragging.kind !== kind || dragging.container !== containerKey || dragging.id === id) return;
      const ids = idsInOrder.slice();
      const from = ids.indexOf(dragging.id), to = ids.indexOf(id);
      if (from === -1 || to === -1) return;
      ids.splice(to, 0, ids.splice(from, 1)[0]);
      commitOrder(kind, idsInOrder, ids);
    });
  };

  const renderSetItem = (en, level, containerKey, idsInOrder, parentName) => {
    const s = en.s;
    const item = el("div", "set-item" + (level === 1 ? " child" : level === 2 ? " variant" : "") + (state.activeType === "set" && s.set_id === state.activeId ? " active" : ""));
    const pct = s.summary?.validated_pct ?? 0;
    const oxHtml = s.standard?.ok ? `<span class="ox ox-o">O</span>` : `<span class="ox ox-x">X</span>`;
    const checked = state.selectedForGroup.has(s.set_id);
    const nVar = en.variants ? en.variants.length : 0;
    const open = searching || !state.treeClosedVariants.has("v:" + s.set_id);
    const dup = s.duplicate_of
      ? `<span class="dup-badge" title="F/R/Probe 서열이 ${escHtml(s.duplicate_of)}와 같습니다 (이름만 다른 중복 등록)">= ${escHtml(s.duplicate_of)}</span>`
      : (s.duplicate_ids && s.duplicate_ids.length ? `<span class="dup-badge" title="같은 서열의 다른 이름: ${escHtml(s.duplicate_ids.join(", "))}">동일 서열 +${s.duplicate_ids.length}</span>` : "");
    const varBadge = nVar ? `<span class="variant-badge" data-toggle="${escHtml(s.set_id)}" title="이 세트 아래의 하위 세트(변이) ${nVar}개 — 눌러서 접기/펼치기">${open ? "▾" : "▸"} 하위 ${nVar}</span>` : "";
    const nm = s.display_name || s.set_id;
    const shortName = level === 2 && parentName && nm.startsWith(parentName) ? "↳ " + nm.slice(parentName.length) + ` <span class="note">${escHtml(nm)}</span>` : escHtml(nm);
    item.innerHTML = `
      <div class="row1">
        <span class="name">${state.groupMode ? `<input type="checkbox" class="set-item-checkbox" data-set-id="${s.set_id}" ${checked ? "checked" : ""}>` : ""}&#8942; ${shortName}</span>
        ${state.groupMode ? "" : `<button class="rename-btn" title="이름 변경" data-set-id="${s.set_id}">&#9998;</button>`}
        <span class="pill pill-${s.segment}">${s.segment}</span>
      </div>
      <div class="stats">validated ${fmtPct(pct)} &middot; ${oxHtml}${s.fluorophore ? ` &middot; <span class="fluor-tag">${escHtml(s.fluorophore)}</span>` : ""} ${varBadge} ${dup}</div>
      <div class="bar"><span style="width:${Math.max(0, Math.min(100, pct))}%"></span></div>
    `;
    if (state.groupMode) {
      item.querySelector(".set-item-checkbox").addEventListener("click", (e) => { e.stopPropagation(); toggleGroupSelection(s.set_id, s.segment); });
      item.addEventListener("click", () => toggleGroupSelection(s.set_id, s.segment));
    } else {
      item.addEventListener("click", () => selectSet(s.set_id));
      item.querySelector(".rename-btn").addEventListener("click", (e) => { e.stopPropagation(); renameSet(s.set_id, nm); });
    }
    const vb = item.querySelector(".variant-badge");
    if (vb) vb.addEventListener("click", (e) => {
      e.stopPropagation();
      const k = "v:" + s.set_id;
      if (state.treeClosedVariants.has(k)) state.treeClosedVariants.delete(k); else state.treeClosedVariants.add(k);
      saveTreeState(); renderSidebar();
    });
    if (level < 2) makeDraggable(item, s.set_id, "set", containerKey, idsInOrder);
    const frag = el("div", null);
    frag.appendChild(item);
    if (nVar && open) {
      const vwrap = el("div", "variant-list");
      for (const v of en.shownVariants) vwrap.appendChild(renderSetItem({ s: v, variants: [], shownVariants: [] }, 2, "", [], nm));
      frag.appendChild(vwrap);
    }
    return frag;
  };

  for (const seg of ["S", "L", ...[...new Set(state.sets.map((x) => x.segment))].filter((x) => x !== "S" && x !== "L")]) {
    const tree = buildSetTree(seg, q);
    if (!tree.length) continue;
    list.appendChild(el("div", null, `<div style="font-size:11px;color:var(--muted);font-weight:700;margin:8px 4px 2px;">SEGMENT ${seg}</div>`));
    const groupIds = tree.filter((t) => t.type === "group").map((t) => t.g.group_id);
    for (const node of tree) {
      if (node.type === "group") {
        const g = node.g, key = "g:" + g.group_id, closed = !searching && state.treeClosed.has(key);
        const item = el("div", "group-item" + (state.activeType === "group" && g.group_id === state.activeId ? " active" : ""));
        item.innerHTML = `
          <div class="row1">
            <span class="name"><span class="tree-chevron" title="접기/펼치기">${closed ? "▸" : "▾"}</span> &#9633; ${escHtml(g.label)} <span style="color:var(--muted);font-weight:400;">(${g.member_count}개 세트 ALL)</span></span>
            <button class="rename-btn" title="그룹 이름 변경" data-group-id="${g.group_id}">&#9998;</button>
            <span class="pill pill-${g.segment}">${g.segment}</span>
          </div>
          <div class="stats">합집합 ${fmtPct(g.summary.validated_pct)} (${g.summary.validated}/${g.summary.total_genomes}) · 교차 포함 <b>${fmtPct(g.summary.validated_with_cross_pct ?? g.summary.validated_pct)}</b></div>
          <div class="bar" title="진한 막대 = 교차 결합 포함, 안쪽 막대 = 세트별 합집합"><span style="width:${Math.max(0, Math.min(100, g.summary.validated_with_cross_pct ?? g.summary.validated_pct))}%"></span></div>
          ${node.childGroups.length ? `<div class="note" style="margin-top:3px;">세트는 아래 그룹에 있습니다: ${node.childGroups.map(escHtml).join(" · ")}</div>` : ""}
        `;
        item.addEventListener("click", () => selectGroup(g.group_id));
        item.querySelector(".tree-chevron").addEventListener("click", (e) => {
          e.stopPropagation();
          if (state.treeClosed.has(key)) state.treeClosed.delete(key); else state.treeClosed.add(key);
          saveTreeState(); renderSidebar();
        });
        item.querySelector(".rename-btn").addEventListener("click", (e) => { e.stopPropagation(); renameGroup(g.group_id, g.label); });
        makeDraggable(item, g.group_id, "group", "groups|" + seg, groupIds);
        list.appendChild(item);
        if (!closed && node.kids.length) {
          const wrap = el("div", "tree-children");
          const ids = node.kids.map((k) => k.s.set_id);
          for (const k of node.kids) wrap.appendChild(renderSetItem(k, 1, "g|" + g.group_id, ids));
          list.appendChild(wrap);
        }
      } else {
        list.appendChild(el("div", null, `<div class="tree-vendor">그룹 없음 · ${escHtml(node.name)} (${node.sets.length})</div>`));
        const wrap = el("div", "source-group");
        const ids = node.sets.map((k) => k.s.set_id);
        for (const k of node.sets) wrap.appendChild(renderSetItem(k, 0, "v|" + seg + "|" + node.name, ids));
        list.appendChild(wrap);
      }
    }
  }
}

function toggleGroupSelection(setId, segment) {
  if (state.selectedForGroup.has(setId)) {
    state.selectedForGroup.delete(setId);
  } else {
    // enforce same-segment selection
    const otherSeg = [...state.selectedForGroup]
      .map((id) => state.sets.find((s) => s.set_id === id)?.segment)
      .find((sg) => sg && sg !== segment);
    if (otherSeg) {
      uiAlert("같은 segment(S 또는 L)의 세트끼리만 묶을 수 있습니다.");
      return;
    }
    state.selectedForGroup.add(setId);
  }
  updateGroupBar();
  renderSidebar();
}

function updateGroupBar() {
  const n = state.selectedForGroup.size;
  $("#groupSelCount").textContent = `${n}개 선택`;
  $("#groupCreateBtn").disabled = n < 2 || !$("#groupLabelInput").value.trim();
}

async function loadSets() {
  state.sets = await fetchJSON("/api/sets");
  renderSidebar();
}

async function loadUiState() {
  state.uiState = await fetchJSON("/api/ui-state");
}

async function saveOrder(kind, order) {
  await fetch("/api/ui-state/order", {
    method: "POST", headers: { "Content-Type": "application/json" },
    body: JSON.stringify({ kind, order }),
  });
}

async function loadGroups() {
  state.groups = await fetchJSON("/api/groups");
  renderSidebar();
}

// ---------------- detail view ----------------

function mmClass(v) {
  if (v === null || v === undefined || v === "") return "mm-na";
  const n = Number(v);
  if (n === 0) return "mm-0";
  if (n === 1) return "mm-1";
  return "mm-2plus";
}
function mmText(v) { return (v === null || v === undefined || v === "") ? "—" : v; }

const ID_SEP = "::"; // matches scripts/group_cross_binding.py's ID_SEP -- "<set_id>::F"/"::R"/"::P"

const FLUOROPHORE_PRESETS = ["FAM", "HEX", "VIC", "CY5", "CY5.5", "ROX", "TEXAS RED", "SUN", "TAMRA"];

function renderFluorophoreOptions(current) {
  const isCustom = current && !FLUOROPHORE_PRESETS.includes(current);
  return `<option value="">(미지정)</option>`
    + FLUOROPHORE_PRESETS.map((f) => `<option value="${f}" ${current === f ? "selected" : ""}>${f}</option>`).join("")
    + `<option value="__custom__" ${isCustom ? "selected" : ""}>직접 입력...</option>`;
}

function renderCitation(citation) {
  if (!citation || !citation.text) return "";
  const link = citation.url
    ? ` <a href="${escHtml(citation.url)}" target="_blank" rel="noopener noreferrer">${escHtml(citation.url)}</a>`
    : "";
  return `<div class="citation-box" style="margin-top:8px;padding:8px 10px;border:1px dashed var(--border);border-radius:8px;font-size:11.5px;color:var(--muted);">
    <b style="color:var(--fg);">출처 (Reference)</b><br>${escHtml(citation.text)}${link}
  </div>`;
}

const TM_TOOLTIP ="Tm(용융 온도, melting temperature): 이 프라이머/프로브와 결합 상대가 이룬 이중나선의 절반이 풀리는 온도(°C). 실제 PCR annealing 온도(보통 55~60°C 부근)보다 충분히 높아야 반응 내내 결합이 유지됩니다.";
const DG_TOOLTIP = "ΔG(깁스 자유에너지, kcal/mol): 이 결합이 얼마나 자발적/안정적으로 형성되는지 나타내는 값. 0에 가까울수록(약한 음수) 결합이 약해 쉽게 떨어지고, 음수가 클수록(예: -20 이하) 결합이 강하고 안정적입니다. mismatch가 0이어도 서열이 짧거나 GC 함량이 낮으면 ΔG가 약해질 수 있어, mismatch 개수만으로는 안 보이는 '결합 안정성'을 보여줍니다.";

const IUPAC_RC = {
  A: "T", C: "G", G: "C", T: "A",
  R: "Y", Y: "R", S: "S", W: "W", K: "M", M: "K",
  B: "V", V: "B", D: "H", H: "D", N: "N",
};
function reverseComplementJS(seq) {
  return seq.toUpperCase().split("").reverse().map((c) => IUPAC_RC[c] || "N").join("");
}
function countIupacMismatches(query, target) {
  let n = 0;
  for (let i = 0; i < Math.min(query.length, target.length); i++) {
    if (!iupacMatches(query[i], target[i])) n++;
  }
  return n;
}

function verdictFor(g) {
  const bad = [];
  if (g.f_mismatch === null || g.f_mismatch === undefined) bad.push("Forward 매치 안됨");
  else if (g.f_mismatch > 1) bad.push(`Forward mismatch ${g.f_mismatch}개`);
  else if (g.f_3p_ok === false) bad.push("Forward 3' 말단에 mismatch 있음");
  if (g.p_mismatch === null || g.p_mismatch === undefined) bad.push("Probe 매치 안됨");
  else if (g.p_mismatch > 0) bad.push(`Probe mismatch ${g.p_mismatch}개(완전일치 필요)`);
  if (g.r_mismatch === null || g.r_mismatch === undefined) bad.push("Reverse 매치 안됨");
  else if (g.r_mismatch > 1) bad.push(`Reverse mismatch ${g.r_mismatch}개`);
  else if (g.r_3p_ok === false) bad.push("Reverse 3' 말단에 mismatch 있음");
  if (bad.length === 0) {
    return { ok: true, text: "✅ 검증 통과 — F/R mismatch ≤1(3' 말단 완전일치), Probe 완전일치 조건을 모두 만족합니다." };
  }
  return { ok: false, text: `⚠️ 검증 실패 — ${bad.join(", ")}` };
}

// Standard-material state: "ok" (validated), "partial" (not validated, but something really binds: F/R paired into an amplicon,
// or any oligo with at most PARTIAL_MM mismatches -- a probe sitting 6 mismatches away is chance, not binding) or "none".
const PARTIAL_MM = 2;
function standardState(g) {
  if (verdictFor(g).ok) return "ok";
  const near = [g.f_mismatch, g.p_mismatch, g.r_mismatch].some((m) => m !== null && m !== undefined && m <= PARTIAL_MM);
  return (g.amplicon_seq || near) ? "partial" : "none";
}

// "Near" binding: not validated (a primer is over the pass tolerance) but it really binds -- probe perfect, F and R both bound, no
// mismatch at either 3' end, and no more than PARTIAL_MM mismatches per primer. Shown as a caution (yellow), not as "unbound" (red).
function nearBinding(g) {
  if (verdictFor(g).ok) return false;
  const has = (m) => m !== null && m !== undefined && m !== "";
  if (!has(g.f_mismatch) || !has(g.r_mismatch) || !has(g.p_mismatch)) return false;
  return Number(g.p_mismatch) === 0 && g.f_3p_ok !== false && g.r_3p_ok !== false &&
    Number(g.f_mismatch) <= PARTIAL_MM && Number(g.r_mismatch) <= PARTIAL_MM;
}
function verdictClass(g, verdict) { return verdict.ok ? "ok" : (nearBinding(g) ? "warn" : "bad"); }

// Server-computed overall-design verdict (webapp/server.py's build_assessment() call,
// same "ok"/"warn"/"bad" classes as narrateSetAssessment() below) rendered as a small pill
// -- e.g. for the all-sets overview table where a full assessment box per row is too heavy.
function evaluationPillHtml(evaluation) {
  if (!evaluation) return "";
  return `<span class="pill pill-${evaluation.cls}">${escHtml(evaluation.label)}</span>`;
}

function assessBindingQuality(g) {
  // Rough, generic TaqMan qPCR design heuristics (not a hard pass/fail) -- flags values
  // far outside typical ranges so a mismatch=0 "perfect" match that would still perform
  // poorly in real qPCR (e.g. a huge amplicon, or a barely-stable primer) doesn't read as
  // unconditionally fine just because the sequence matches.
  const notes = [];
  const hasAmpRange = g.amplicon_start !== null && g.amplicon_start !== undefined && g.amplicon_start !== "" &&
    g.amplicon_end !== null && g.amplicon_end !== undefined && g.amplicon_end !== "";
  if (hasAmpRange) {
    const ampLen = Number(g.amplicon_end) - Number(g.amplicon_start) + 1;
    if (!isNaN(ampLen)) {
      if (ampLen > 300) notes.push({ ok: false, text: `Amplicon 크기 ${ampLen}bp — TaqMan qPCR 권장 범위(약 60~150bp)보다 매우 커서 증폭 효율 저하나 실패 위험이 있습니다.` });
      else if (ampLen > 150) notes.push({ ok: null, text: `Amplicon 크기 ${ampLen}bp — qPCR 권장 범위(약 60~150bp)보다 다소 크지만 동작 가능한 범위입니다.` });
      else notes.push({ ok: true, text: `Amplicon 크기 ${ampLen}bp — qPCR에 적합한 범위입니다.` });
    }
  }
  const checkOligo = (label, tm, dg, isProbe) => {
    if (tm !== null && tm !== undefined && tm !== "") {
      const tmNum = Number(tm);
      const loTm = isProbe ? 63 : 55, hiTm = isProbe ? 72 : 65;
      if (tmNum < loTm) notes.push({ ok: false, text: `${label} Tm ${tmNum.toFixed(1)}°C — 권장 범위(${loTm}~${hiTm}°C)보다 낮아 표준 annealing 온도에서 결합이 불안정할 수 있습니다.` });
      else if (tmNum > hiTm) notes.push({ ok: null, text: `${label} Tm ${tmNum.toFixed(1)}°C — 권장 범위(${loTm}~${hiTm}°C)보다 높은 편입니다.` });
      else notes.push({ ok: true, text: `${label} Tm ${tmNum.toFixed(1)}°C — 정상 범위입니다.` });
    }
    if (dg !== null && dg !== undefined && dg !== "") {
      const dgNum = Number(dg);
      if (dgNum > -14) notes.push({ ok: false, text: `${label} ΔG ${dgNum.toFixed(2)} kcal/mol — 결합 자유에너지가 약해(덜 음수) 결합이 불안정하거나 비특이적일 위험이 있습니다.` });
    }
  };
  checkOligo("Forward", g.f_tm, g.f_dg, false);
  checkOligo("Probe", g.p_tm, g.p_dg, true);
  checkOligo("Reverse", g.r_tm, g.r_dg, false);
  return notes;
}

function renderBindingQuality(g) {
  const notes = assessBindingQuality(g);
  if (!notes.length) return "";
  const icon = (ok) => (ok === true ? "✅" : ok === false ? "⚠️" : "ℹ️");
  const cls = (ok) => (ok === true ? "ok" : ok === false ? "bad" : "warn");
  return `
    <div class="amp-quality">
      <div class="amp-quality-title">결합 품질 평가 (참고용 — amplicon 크기/Tm/ΔG 기반)</div>
      ${notes.map((n) => `<div class="amp-quality-row ${cls(n.ok)}">${icon(n.ok)} ${n.text}</div>`).join("")}
    </div>`;
}

function highlightAmplicon(g, setData, labels) {
  labels = labels || { f: "FORWARD", p: "PROBE", r: "REVERSE" };
  let seq = g.amplicon_seq || "";
  if (!seq) {
    // F/R never paired up, so there's no amplicon to align against -- but the probe (searched
    // independently) or a partial F/R signal may still exist, and exactly how far off a
    // mismatch is is the point of this view, so show whatever numbers are available instead
    // of just "no sequence".
    const noAmpRow = (roleClass, label, mm, tm, dg) => `
      <div class="amp-role-row">
        <span class="amp-role-tag ${roleClass}">${label}</span>
        <span class="msa-mm-tag ${mmClass(mm)}">mismatch ${mmText(mm)}</span>
        ${tm !== null && tm !== undefined && tm !== "" ? `<span class="amp-thermo mono" title="${TM_TOOLTIP}">Tm ${Number(tm).toFixed(1)}°C</span>` : ""}
        ${dg !== null && dg !== undefined && dg !== "" ? `<span class="amp-thermo mono" title="${DG_TOOLTIP}">ΔG ${Number(dg).toFixed(2)} kcal/mol</span>` : ""}
      </div>`;
    const verdict = verdictFor(g);
    return `
      <div class="amp-verdict ${verdictClass(g, verdict)}">${verdict.text}${nearBinding(g) ? " <span class=\"note\">(결합은 확인됨 — 3′ 끝 정상, 프라이머 mismatch가 검증 기준만 넘음)</span>" : ""}</div>
      <div class="amp-role-summary">
        ${noAmpRow("f", labels.f, g.f_mismatch, g.f_tm, g.f_dg)}
        ${noAmpRow("p", labels.p, g.p_mismatch, g.p_tm, g.p_dg)}
        ${noAmpRow("r", labels.r, g.r_mismatch, g.r_tm, g.r_dg)}
      </div>
      <p class="note" style="margin-top:6px;">F/R이 짝을 이루지 못해 amplicon 서열을 얻지 못했습니다 (정렬 표시는 amplicon이 있어야 가능합니다).</p>
      ${renderBindingQuality(g)}`;
  }
  const rawSeq = seq;
  const n = rawSeq.length;
  const ampStart = Number(g.amplicon_start);
  const ampEnd = Number(g.amplicon_end);
  const pStartAbs = g.probe_start !== null && g.probe_start !== undefined && g.probe_start !== "" ? Number(g.probe_start) : null;
  const pSeq = setData?.P_seq || "";
  const fSeq = setData?.F_seq || "";
  const rSeq = setData?.R_seq || "";
  const fLen = fSeq.length;
  const rLen = rSeq.length;

  // mfeprimer's amplicon-extraction convention assumes the genome record is deposited on
  // the same strand as the reference used to design F/R, so F always sits at the
  // amplicon's own start and R (reverse-complemented) at its end -- true for most genomes,
  // but some accessions are deposited under the opposite strand convention, which flips
  // this (confirmed for OQ919513/VendorD_GPC_PREVIOUS1_LASV1_P2: R matched the amplicon's own
  // first bases directly and RC(F) matched its last bases -- exactly backwards from the
  // assumed convention, which showed up as a garbled/misaligned view for a genome the
  // mismatch table itself correctly reported as a perfect 0/0/0 match). Self-detect which
  // orientation actually fits (same idea as the probe strand check below) and normalize
  // the *displayed* sequence to it, rather than assuming the convention always holds.
  let flipped = false;
  if (fLen && rLen && n) {
    const rRc = reverseComplementJS(rSeq);
    const directMM = countIupacMismatches(fSeq, rawSeq.slice(0, fLen)) + countIupacMismatches(rRc, rawSeq.slice(Math.max(0, n - rLen)));
    const flippedCheck = reverseComplementJS(rawSeq);
    const flippedMM = countIupacMismatches(fSeq, flippedCheck.slice(0, fLen)) + countIupacMismatches(rRc, flippedCheck.slice(Math.max(0, n - rLen)));
    flipped = flippedMM < directMM;
  }
  seq = flipped ? reverseComplementJS(rawSeq) : rawSeq;

  // F/R coordinate LABELS always describe true genome positions -- when the display got
  // flipped, F is actually the oligo sitting at the genome's amp_end side, and vice versa.
  const fCoordStart = flipped ? (!isNaN(ampEnd) ? ampEnd - fLen + 1 : null) : (!isNaN(ampStart) ? ampStart : null);
  const fCoordEnd = flipped ? (!isNaN(ampEnd) ? ampEnd : null) : (!isNaN(ampStart) ? ampStart + fLen - 1 : null);
  const rCoordStart = flipped ? (!isNaN(ampStart) ? ampStart : null) : (!isNaN(ampEnd) ? ampEnd - rLen + 1 : null);
  const rCoordEnd = flipped ? (!isNaN(ampStart) ? ampStart + rLen - 1 : null) : (!isNaN(ampEnd) ? ampEnd : null);

  // The probe's own absolute genome coordinate (used for its label) is unaffected by the
  // display flip -- only its LOCAL character offset within the now-possibly-reoriented
  // `seq` string needs to account for it, mirrored around the amplicon span (reversing an
  // n-length array maps old index i -> n-1-i, so a [a,b) span becomes [n-b, n-a)).
  let pLocalStart = null, pLocalEnd = null;
  if (pStartAbs !== null && !isNaN(ampStart) && !isNaN(ampEnd)) {
    if (flipped) {
      pLocalEnd = ampEnd - pStartAbs + 1;
      pLocalStart = pLocalEnd - pSeq.length;
    } else {
      pLocalStart = pStartAbs - ampStart;
      pLocalEnd = pLocalStart + pSeq.length;
    }
  }

  const rStart = Math.max(0, n - rLen);
  const role = new Array(n).fill(null);
  for (let k = 0; k < Math.min(fLen, n); k++) role[k] = "f";
  for (let k = rStart; k < n; k++) role[k] = "r";

  let pStrand = "+";
  const hasProbe = pLocalStart !== null && pLocalEnd !== null && pLocalStart >= 0 && pLocalEnd <= n && pLocalStart < pLocalEnd;
  if (hasProbe) {
    const pHit = seq.slice(pLocalStart, pLocalEnd);
    const rc = countIupacMismatches(pSeq, reverseComplementJS(pHit));
    if (rc < countIupacMismatches(pSeq, pHit)) pStrand = "-";
    for (let k = pLocalStart; k < pLocalEnd; k++) role[k] = "p";
  }

  // One continuous, char-for-char aligned 3-line block (fixed-width font is mandatory here
  // -- see .amp-align in style.css) spanning the *whole* amplicon: the primer/probe's own
  // sequence sitting exactly where it binds, a ":" wherever that base satisfies the
  // primer's IUPAC code, and the amplicon itself below -- so binding is visible in context
  // instead of three disconnected boxes.
  const rDisplay = reverseComplementJS(rSeq);
  const pDisplay = pStrand === "-" ? reverseComplementJS(pSeq) : pSeq;
  const primerChars = new Array(n).fill("-");
  const colonChars = new Array(n).fill(" ");
  const dirChars = new Array(n).fill(" ");
  for (let k = 0; k < n; k++) {
    if (role[k] === "f") primerChars[k] = fSeq[k];
    else if (role[k] === "r") primerChars[k] = rDisplay[k - rStart];
    else if (role[k] === "p") primerChars[k] = pDisplay[k - pLocalStart];
    if (role[k]) colonChars[k] = iupacMatches(primerChars[k], seq[k]) ? ":" : "X";
  }
  // Mismatch colours: RED = any probe mismatch, or a primer mismatch that matters (within the primer's 3' end window --
  // config three_prime_consecutive_match -- or more mismatches than the verdict tolerates); YELLOW = a primer mismatch
  // away from the 3' end that does not break binding.
  const n3 = Number(state.params?.three_prime_consecutive_match) || 2;
  const spanOf = { f: [0, Math.min(fLen, n)], r: [rStart, n] };
  const mmSeverity = new Array(n).fill(null);
  for (const r of ["f", "r"]) {
    const [a, b] = spanOf[r];
    const bad = [];
    for (let k = a; k < b; k++) if (colonChars[k] === "X") bad.push(k);
    const tooMany = bad.length > PARTIAL_MM;
    for (const k of bad) {
      const in3p = r === "f" ? k >= b - n3 : k < a + n3;               // F's 3' end is its right end, R (shown reverse-complemented) its left end
      mmSeverity[k] = (in3p || tooMany) ? "red" : "yellow";
    }
  }
  for (let k = 0; k < n; k++) if (role[k] === "p" && colonChars[k] === "X") mmSeverity[k] = "red";
  const cell = (k, ch) => (mmSeverity[k] ? `<span class="mm-${mmSeverity[k]}">${ch}</span>` : ch);
  // 5'/3' direction markers over each span, drawn in the amplicon's own left-to-right
  // display order -- a "-" strand oligo (R, and Probe when it binds the other strand) is
  // shown reverse-complemented so its letters line up with the amplicon above, so its own
  // 5' end visually sits on the RIGHT (arrows point left) even though the span itself
  // still reads left-to-right on screen.
  const writeDirection = (start, end, strand) => {
    const width = end - start;
    if (width < 4) return; // too narrow for "5'"..arrows.."3'" to mean anything
    const arrow = strand === "-" ? "<" : ">";
    for (let k = start; k < end; k++) dirChars[k] = arrow;
    const left = strand === "-" ? "3'" : "5'";
    const right = strand === "-" ? "5'" : "3'";
    dirChars[start] = left[0]; dirChars[start + 1] = left[1];
    dirChars[end - 2] = right[0]; dirChars[end - 1] = right[1];
  };
  writeDirection(0, Math.min(fLen, n), "+");
  writeDirection(rStart, n, "-");
  if (hasProbe) writeDirection(pLocalStart, pLocalEnd, pStrand);

  const wrapByRole = (chars) => {
    let out = "", i = 0;
    while (i < n) {
      let j = i;
      while (j < n && role[j] === role[i]) j++;
      const chunk = chars.slice(i, j).map((c, off) => cell(i + off, c)).join("");
      out += role[i] ? `<mark class="${role[i]}">${chunk}</mark>` : chunk;
      i = j;
    }
    return out;
  };

  const roleRow = (roleClass, label, strand, mm, start, end, tm, dg) => `
    <span class="amp-role-tag ${roleClass}">${label}${strand ? ` (${strandLabel(strand)})` : ""}</span>
    ${start != null && !isNaN(start) ? `<span class="mono">${start}-${end}</span>` : "<span class='note'>매치 안됨</span>"}
    <span class="msa-mm-tag ${mmClass(mm)}">mismatch ${mmText(mm)}</span>
    ${tm !== null && tm !== undefined && tm !== "" ? `<span class="amp-thermo mono" title="${TM_TOOLTIP}">Tm ${Number(tm).toFixed(1)}°C</span>` : ""}
    ${dg !== null && dg !== undefined && dg !== "" ? `<span class="amp-thermo mono" title="${DG_TOOLTIP}">ΔG ${Number(dg).toFixed(2)} kcal/mol</span>` : ""}`;

  const verdict = verdictFor(g);
  return `
    <div class="amp-verdict ${verdictClass(g, verdict)}">${verdict.text}${nearBinding(g) ? " <span class=\"note\">(결합은 확인됨 — 3′ 끝 정상, 프라이머 mismatch가 검증 기준만 넘음)</span>" : ""}</div>
    <div class="amp-role-summary">
      <div class="amp-role-row">${roleRow("f", labels.f, flipped ? "-" : "+", g.f_mismatch, fCoordStart, fCoordEnd, g.f_tm, g.f_dg)}</div>
      <div class="amp-role-row">${roleRow("p", labels.p, hasProbe ? pStrand : null, g.p_mismatch, pStartAbs, hasProbe ? pStartAbs + pSeq.length - 1 : null, g.p_tm, g.p_dg)}</div>
      <div class="amp-role-row">${roleRow("r", labels.r, flipped ? "+" : "-", g.r_mismatch, rCoordStart, rCoordEnd, g.r_tm, g.r_dg)}</div>
    </div>
    <div class="amp-align">
      <div class="amp-align-row amp-align-dir">${dirChars.join("")}</div>
      <div class="amp-align-row">${wrapByRole(primerChars)}</div>
      <div class="amp-align-row amp-align-colon">${colonChars.map((c, k) => cell(k, c)).join("")}</div>
      <div class="amp-align-row">${wrapByRole(seq.split(""))}</div>
    </div>
    ${mmSeverity.some(Boolean) ? `<div class="amp-mm-legend note"><span class="mm-red">빨강</span> 프로브 mismatch, 또는 프라이머 3′ 말단(끝 ${n3}염기) mismatch · ${PARTIAL_MM}개 초과 &nbsp; <span class="mm-yellow">노랑</span> 프라이머 mismatch (3′ 말단 밖, ${PARTIAL_MM}개 이하 — 결합은 됨)</div>` : ""}
    ${renderBindingQuality(g)}`;
}

function standardSummaryTable(standardGenomes) {
  return `<table class="mini-reference-table"><tr><th>표준물질</th><th>판정</th><th>Forward</th><th>Probe</th><th>Reverse</th><th>amplicon</th></tr>${
    standardGenomes.map((g) => {
      const v = verdictFor(g);
      const amp = g.amplicon_start != null && g.amplicon_start !== "" ? `${g.amplicon_start}-${g.amplicon_end}` : "&mdash;";
      return `<tr><td class="mono">${escHtml(g.genome_id)}</td>
        <td><span class="pill ${v.ok ? "pill-ok" : standardState(g) === "partial" ? "pill-warn" : "pill-bad"}">${v.ok ? "O 매치" : standardState(g) === "partial" ? "△ 일부 매치" : "X 불일치"}</span></td>
        <td>${mmText(g.f_mismatch)}${g.f_3p_ok === false ? " (3′ 불일치)" : ""}</td><td>${mmText(g.p_mismatch)}</td>
        <td>${mmText(g.r_mismatch)}${g.r_3p_ok === false ? " (3′ 불일치)" : ""}</td><td>${amp}</td></tr>`;
    }).join("")}</table>`;
}

// Pop-up for one set: every standard material with its verdict, mismatch counts and the binding alignment.
function openStandardModal(setData) {
  const sg = setData.standard_genomes || [];
  const ok = sg.filter((g) => verdictFor(g).ok).length;
  const part = sg.filter((g) => standardState(g) === "partial").length;
  return openDetailModal({
    title: `표준물질(KDCA) 매치 &mdash; ${escHtml(setData.display_name || setData.set_id)}`,
    subtitle: sg.length ? `${sg.length}개 중 ${ok}개 일치${part ? ` · ${part}개 일부 매치` : ""} · 판정: F/R mismatch&le;1(3′ 끝 일치), Probe 완전일치` : "",
    html: sg.length ? standardSummaryTable(sg) + renderStandardDetail(sg, setData, { skipUnmatched: true })
      : "<p class='note'>이 세트의 segment에는 등록된 표준물질이 없습니다 (표준물질은 선택 사항입니다).</p>",
  });
}

// Pop-up for the overview tile: every set x every standard material.
async function openStandardOverviewModal() {
  const ids = (state.sets || []).map((s) => s.set_id);
  const details = (await Promise.all(ids.map((id) => fetchJSON(`/api/sets/${encodeURIComponent(id)}`).catch(() => null)))).filter(Boolean);
  const names = [...new Set(details.flatMap((d) => (d.standard_genomes || []).map((g) => g.genome_id)))];
  if (!names.length) {
    return openDetailModal({ title: "표준물질 매치", subtitle: "", html: "<p class='note'>등록된 표준물질이 없습니다 (표준물질은 선택 사항입니다).</p>" });
  }
  const head = `<tr><th>프라이머 세트</th>${names.map((n) => `<th class="mono">${escHtml(n)}</th>`).join("")}<th>일치</th></tr>`;
  const body = details.map((d) => {
    const by = Object.fromEntries((d.standard_genomes || []).map((g) => [g.genome_id, g]));
    let okCount = 0, total = 0;
    const cells = names.map((n) => {
      const g = by[n];
      if (!g) return "<td>&mdash;</td>";
      const v = verdictFor(g), st = standardState(g); total += 1; okCount += v.ok ? 1 : 0;
      return `<td><span class="pill ${st === "ok" ? "pill-ok" : st === "partial" ? "pill-warn" : "pill-bad"}">${st === "ok" ? "O" : st === "partial" ? "△" : "X"}</span> <span class="note">${g.f_mismatch == null && g.r_mismatch == null ? "결합 없음" : `F${g.f_mismatch ?? "?"} P${g.p_mismatch ?? "?"} R${g.r_mismatch ?? "?"}`}</span></td>`;
    }).join("");
    return `<tr><td>${escHtml(d.display_name || d.set_id)} <span class="note">${escHtml(d.set_id)}</span></td>${cells}<td>${total ? `${okCount} / ${total}` : "&mdash;"}</td></tr>`;
  }).join("");
  return openDetailModal({
    title: "표준물질 매치 &mdash; 전체 세트",
    subtitle: "O = F/R mismatch&le;1(3′ 끝 일치) 그리고 Probe 완전일치 · 한 세트의 자세한 결합 정렬은 그 세트의 \"표준물질(KDCA) 매치\" 카드를 누르세요",
    html: `<div style="overflow:auto;"><table class="mini-reference-table">${head}${body}</table></div>`,
  });
}

document.addEventListener("click", async (e) => {
  const run = e.target.closest("#analyzeSetBtn, #reanalyzeSetBtn");
  if (run && state.detailData) {
    const force = run.id === "reanalyzeSetBtn", id = state.detailData.set_id;
    if (force && !(await uiConfirm(`세트 ${id}의 결과를 처음부터 다시 계산합니다. 계속할까요?`, { title: "이 세트 재분석", okLabel: "재분석" }))) return;
    startPipelineJob(`/api/sets/${encodeURIComponent(id)}/analyze${force ? "?force=1" : ""}`, force ? `세트 ${id} 재분석` : `세트 ${id} 분석`);
    return;
  }
  const grun = e.target.closest("#analyzeGroupBtn, #reanalyzeGroupBtn");
  if (grun && state.detailGroup) {
    const force = grun.id === "reanalyzeGroupBtn", id = state.detailGroup.group_id;
    if (force && !(await uiConfirm(`그룹 ${id}의 조합 결과(합집합·교차 결합·표준물질/참조유전체 조합)를 처음부터 다시 계산합니다. 계속할까요?`, { title: "이 그룹 재분석", okLabel: "재분석" }))) return;
    startPipelineJob(`/api/groups/${encodeURIComponent(id)}/analyze${force ? "?force=1" : ""}`, force ? `그룹 ${id} 재분석` : `그룹 ${id} 분석`);
    return;
  }
  if (e.target.closest("#stdTile")) openStandardOverviewModal();
  else if (e.target.closest("#stdMatchCard") && state.detailData) openStandardModal(state.detailData);
  else if (e.target.closest("#bindingCard") && state.detailData) openBindingModal(state.detailData);
});
document.addEventListener("keydown", (e) => {
  if (e.key !== "Enter" || !document.activeElement) return;
  if (document.activeElement.id === "stdTile") openStandardOverviewModal();
  else if (document.activeElement.id === "stdMatchCard" && state.detailData) openStandardModal(state.detailData);
  else if (document.activeElement.id === "bindingCard" && state.detailData) openBindingModal(state.detailData);
});

// Pop-up for the "참조유전체 결합위치" card: where F / Probe / R bind on the reference genome, in detail.
function bindingTrackSvg(binding, gene) {
  const hits = (binding || []).filter((h) => h.found);
  if (!hits.length) return "";
  const pos = hits.flatMap((h) => [Number(h.start), Number(h.end)]);
  const lo = Math.min(...pos), hi = Math.max(...pos), pad = Math.max(10, Math.round((hi - lo) * 0.08));
  const a = lo - pad, b = hi + pad, W = 700, x = (p) => 10 + ((p - a) / (b - a)) * (W - 20);
  const color = { F: "#4C78A8", P: "#F2A65A", R: "#54A24B" };
  const rows = hits.map((h, i) => {
    const y = 40 + i * 24, x1 = x(Number(h.start)), x2 = Math.max(x(Number(h.end)), x1 + 3);
    return `<rect x="${x1.toFixed(1)}" y="${y}" width="${(x2 - x1).toFixed(1)}" height="14" rx="3" fill="${color[h.role] || "#888"}"></rect>
      ${x2 > W - 150 ? `<text x="${(x1 - 6).toFixed(1)}" y="${y + 11}" font-size="11" text-anchor="end" fill="currentColor">` : `<text x="${(x2 + 6).toFixed(1)}" y="${y + 11}" font-size="11" fill="currentColor">`}${h.role} ${h.start}-${h.end} (${h.strand})</text>`;
  }).join("");
  const ampLo = gene && gene.amplicon_start ? Number(gene.amplicon_start) : null;
  const ampHi = gene && gene.amplicon_end ? Number(gene.amplicon_end) : null;
  const amp = ampLo && ampHi ? `<line x1="${x(ampLo).toFixed(1)}" x2="${x(ampHi).toFixed(1)}" y1="26" y2="26" stroke="currentColor" stroke-width="2"></line>
      <text x="${x(ampLo).toFixed(1)}" y="18" font-size="11" fill="currentColor">amplicon ${ampLo}-${ampHi} (${ampHi - ampLo + 1} bp)</text>` : "";
  return `<svg viewBox="0 0 ${W} ${60 + hits.length * 24}" style="width:100%;max-width:${W}px;color:var(--fg);">${amp}${rows}</svg>`;
}

function openBindingModal(setData) {
  const binding = setData.binding || [];
  const gene = setData.gene || {};
  const ref = (binding.find((h) => h.found) || {}).chrom || gene.chrom || "";
  const roleName = { F: "Forward", P: "Probe", R: "Reverse" };
  const rows = binding.map((h) => h.found
    ? `<tr><td>${roleName[h.role] || h.role}</td><td class="mono">${escHtml(h.oligo_id)}</td><td class="mono">${escHtml(h.seq)}</td>
        <td>${escHtml(h.chrom)}:${h.start}-${h.end}</td><td>${h.strand === "+" ? "정방향(+)" : "역방향(&minus;)"}</td><td>${Math.abs(Number(h.end) - Number(h.start)) + 1} bp</td>
        <td><span class="mm ${mmClass(Number(h.mismatches))}" style="padding:1px 8px;border-radius:5px;">${h.mismatches}</span></td></tr>`
    : `<tr><td>${roleName[h.role] || h.role}</td><td class="mono">${escHtml(h.oligo_id)}</td><td class="mono">${escHtml(h.seq)}</td>
        <td colspan="4" style="color:var(--bad);font-weight:600;">NO MATCH &mdash; 참조유전체에서 결합을 찾지 못했습니다</td></tr>`).join("");
  const refGenome = (setData.genomes || []).find((g) => g.genome_id === ref);
  const geneLine = gene.target_gene && gene.target_gene !== "NA"
    ? `<p>표적 유전자: <b>${escHtml(gene.target_gene)}</b>${gene.amplicon_size ? ` · amplicon ${escHtml(String(gene.amplicon_size))} bp (${gene.amplicon_start}-${gene.amplicon_end})` : ""}</p>` : "";
  const found = binding.filter((h) => h.found).length;
  return openDetailModal({
    title: `참조유전체 결합위치 &mdash; ${escHtml(setData.display_name || setData.set_id)}`,
    subtitle: `참조유전체 ${escHtml(ref || "—")} · ${binding.length}개 중 ${found}개 결합`,
    html: `${geneLine}${bindingTrackSvg(binding, gene)}
      <table class="mini-reference-table"><tr><th>역할</th><th>ID</th><th>서열 (5′→3′)</th><th>위치</th><th>가닥</th><th>길이</th><th>mismatch</th></tr>${rows}</table>
      <p class="note">좌표는 참조유전체 기준 1-based(양 끝 포함)이며, 역방향(&minus;)은 참조유전체의 반대 가닥에 결합한다는 뜻입니다.</p>
      ${refGenome && refGenome.amplicon_seq ? `<h4 style="margin:14px 0 6px;">참조유전체의 결합 정렬</h4>${highlightAmplicon(refGenome, setData)}` : ""}`,
  });
}

function renderStandardDetail(standardGenomes, setData, opts) {
  if (!standardGenomes || !standardGenomes.length) {
    return "<p class='note'>표준물질 매치 데이터가 없습니다 (분석 다시하기가 필요할 수 있습니다).</p>";
  }
  // Only the standards that match are shown with their binding alignment; the others are reduced to one line each
  // (name + why), folded away when something does match.
  const matched = standardGenomes.filter((g) => standardState(g) === "ok");
  const partial = standardGenomes.filter((g) => standardState(g) === "partial");
  const unmatched = standardGenomes.filter((g) => standardState(g) === "none");
  const block = (g, pill) => `
    <div class="standard-material-block">
      <div class="standard-material-head">
        <span class="standard-material-id mono">${escHtml(g.genome_id)}</span>
        ${pill}
      </div>
      ${highlightAmplicon(g, setData)}
    </div>`;
  const matchedHtml = matched.map((g) => block(g, `<span class="pill pill-ok">O 매치</span>`)).join("");
  const mmTxt = (m) => (m === null || m === undefined ? "—" : m);
  const partialHtml = partial.map((g) => block(g, `<span class="pill pill-warn" title="${escHtml(verdictFor(g).text)}">△ 일부 매치</span> <span class="note">F ${mmTxt(g.f_mismatch)} · Probe ${mmTxt(g.p_mismatch)} · R ${mmTxt(g.r_mismatch)} (mismatch)</span>`)).join("");
  const line = (g) => `<li><span class="mono">${escHtml(g.genome_id)}</span> <span class="pill pill-bad">X</span> <span class="note">${escHtml(verdictFor(g).text.replace(/^[^\s]+\s*검증 실패 — /, ""))}</span></li>`;
  const unmatchedHtml = (unmatched.length && !(opts && opts.skipUnmatched))
    ? (matched.length
      ? `<details class="std-unmatched"><summary class="note">결합이 거의 없는 표준물질 ${unmatched.length}개 (펼치기)</summary><ul class="std-unmatched-list">${unmatched.map(line).join("")}</ul></details>`
      : `<ul class="std-unmatched-list">${unmatched.map(line).join("")}</ul>`) : "";
  return (matched.length ? matchedHtml : "<p class='note'>완전히 일치하는 표준물질이 없습니다.</p>")
    + (partial.length ? `<h4 class="std-partial-h">일부만 매치되는 표준물질 ${partial.length}개 <span class="note">(검증 기준은 못 맞췄지만 결합은 확인됨)</span></h4>${partialHtml}` : "") + unmatchedHtml;
}

function renderGenomeTable(genomes, setData) {
  const search = ($("#genomeSearch")?.value || "").toLowerCase();
  let rows = genomes.filter((g) => {
    if (state.problemsOnly && comboRowStatus(g).level === "ok") return false;
    if (search) {
      const hay = `${g.genome_id} ${g.lineage} ${g.country} ${g.region}`.toLowerCase();
      if (!hay.includes(search)) return false;
    }
    return true;
  });
  const { key, dir } = state.genomeSort;
  if (key) {
    rows = rows.slice().sort((a, b) => {
      let av = a[key], bv = b[key];
      if (av === null || av === undefined) av = -1;
      if (bv === null || bv === undefined) bv = -1;
      if (typeof av === "string") return dir * av.localeCompare(bv);
      return dir * (av - bv);
    });
  }
  $("#genomeCount").textContent = `${rows.length} / ${genomes.length}개 게놈`;

  const tbody = $("#genomeTableBody");
  tbody.innerHTML = "";
  const frag = document.createDocumentFragment();
  for (const g of rows.slice(0, 2000)) {
    const status = comboRowStatus(g);
    const tr = el("tr", "genome-row" + (status.level !== "ok" ? " problem-row" : ""));
    tr.dataset.genomeId = g.genome_id;
    tr.innerHTML = `
      <td><span class="combo-status combo-status-${status.level}">${status.icon} ${status.label}</span></td>
      <td class="mono">${g.genome_id}</td>
      <td>${g.lineage}</td>
      <td>${g.country || ""}</td>
      <td class="mm ${mmClass(g.f_mismatch)}">${mmText(g.f_mismatch)}</td>
      <td class="mm ${mmClass(g.p_mismatch)}">${mmText(g.p_mismatch)}</td>
      <td class="mm ${mmClass(g.r_mismatch)}">${mmText(g.r_mismatch)}</td>
    `;
    tr.addEventListener("click", () => openGenomeModal(g, setData));
    frag.appendChild(tr);
  }
  tbody.appendChild(frag);
}

// ---------------- MSA column (rectangular/box) selection, VS Code Alt+drag style ----------------
// A single module-level drag state (rather than per-call closures) plus ONE set of global
// mousemove/mouseup listeners registered once at load time -- wireColumnSelection() is
// re-invoked every time a set's detail view (and its msaCard) is rebuilt, so listeners
// scoped to a specific msaCard/closure would otherwise pile up on `document` forever as
// the user browses between sets.
let msaColSelectDrag = null; // { msaCard, startCol } while a drag is in progress

function msaCharWidth(msaCard) {
  // width of the TEXT (not the padded box) of a full-length row / its number of characters
  const seqEl = msaCard.querySelector(".msa-consensus .msa-seq") || msaCard.querySelector(".msa-ruler .msa-seq");
  if (!seqEl || !seqEl.textContent.length) return 7;
  const rng = document.createRange();
  rng.selectNodeContents(seqEl);
  return rng.getBoundingClientRect().width / seqEl.textContent.length;
}

function msaTextLeft(seqEl) {            // x where the first character starts (the .msa-seq box has left padding)
  return seqEl.getBoundingClientRect().left + (parseFloat(getComputedStyle(seqEl).paddingLeft) || 0);
}

function msaColAt(msaCard, clientX) {
  const seqEl = msaCard.querySelector(".msa-ruler .msa-seq");
  if (!seqEl) return 0;
  const cw = msaCharWidth(msaCard);
  return Math.max(0, Math.floor((clientX - msaTextLeft(seqEl)) / cw));
}

function msaUpdateColOverlay(msaCard, c1, c2) {
  const outer = msaCard.querySelector(".msa-outer");
  const seqEl = msaCard.querySelector(".msa-ruler .msa-seq");
  if (!outer || !seqEl) return;
  let ov = outer.querySelector(".msa-col-select-overlay");
  if (!ov) {
    ov = el("div", "msa-col-select-overlay");
    outer.appendChild(ov);
  }
  const cw = msaCharWidth(msaCard);
  const lo = Math.min(c1, c2), hi = Math.max(c1, c2);
  const outerRect = outer.getBoundingClientRect();
  const left = msaTextLeft(seqEl) - outerRect.left + outer.scrollLeft + lo * cw;
  ov.style.left = `${left}px`;
  ov.style.width = `${Math.max(cw, (hi - lo + 1) * cw)}px`;
  ov.style.height = `${outer.scrollHeight}px`;
  ov.hidden = false;
}

function msaShowColReadout(msaCard, c1, c2) {
  const readout = msaCard.querySelector("#msaColSelectReadout");
  if (!readout) return;
  const lo = Math.min(c1, c2), hi = Math.max(c1, c2);
  readout.hidden = false;
  readout.querySelector(".msa-col-select-range").textContent = `${lo + 1}-${hi + 1} (${hi - lo + 1}bp)`;
  readout.dataset.start = lo;
  readout.dataset.end = hi;
}

function msaClearColSelection(msaCard) {
  const ov = msaCard.querySelector(".msa-col-select-overlay");
  if (ov) ov.hidden = true;
  const readout = msaCard.querySelector("#msaColSelectReadout");
  if (readout) readout.hidden = true;
}

// Probe auto-design (scripts/probe_designer.py): GET shows the last saved run, POST runs it again (about a second).
function renderProbeDesign(r, ctx = { vals: new Map(), running: null, setData: {} }) {
  if (!r || !r.generated) return "<p class='note'>아직 실행하지 않았습니다. 오른쪽 위 \"Probe 설계 실행\"을 누르면 이 세트의 정렬 consensus에서 후보를 만들고 전체 population으로 다시 검증합니다.</p>";
  const reg = r.registered || {};
  const n = r.msa_n || 0, pct = (v) => (n ? ` (${(100 * v / n).toFixed(1)}%)` : "");
  const rows = (r.candidates || []).map((c) => {
    const better = reg.msa_cov0 != null && c.msa_cov0 > reg.msa_cov0;
    const val = ctx.vals.get(c.seq), running = ctx.running && ctx.running.seq === c.seq;
    const valRow = (val || running) ? `<tr class="probe-val-row"><td colspan="13" style="padding:0 0 10px 24px;">${running ? `<p class="note">⏳ 검증 중: ${escHtml(ctx.running.text || "")}</p>` : ""}${val ? probeValidationHtml(val, ctx.setData, true) : ""}</td></tr>` : "";
    return `<tr><td>${c.rank}</td>
      <td class="mono">${escHtml(c.seq)} <button type="button" class="btn btn-secondary probe-copy" data-seq="${escHtml(c.seq)}" title="서열 복사" style="padding:0 6px;">&#128203;</button></td>
      <td>${({ own: "자체 규칙", primer3: "Primer3", "own+primer3": "자체+Primer3" })[c.source] || c.source}${c.p3_penalty != null ? `<div class="note">penalty ${c.p3_penalty} · hairpin ${c.p3_hairpin}&deg;C</div>` : ""}</td><td>${c.strand === "+" ? "정방향 (5′→3′, 센스)" : "역방향 (상보가닥)"}</td><td>${c.length}</td><td>${c.n_iupac ?? 0}</td><td>${c.tm}${c.tm_p3 != null ? ` <span class="note">(${c.tm_p3})</span>` : ""}</td><td>${c.gc}%</td>
      <td>${c.col_start}-${c.col_end}</td>
      <td>${c.msa_cov0}${pct(c.msa_cov0)} ${better ? '<span class="pill pill-ok">등록 Probe보다 높음</span>' : ""}</td><td>${c.msa_cov1}${pct(c.msa_cov1)}</td>
      <td>${c.mfe_hits != null ? `${c.mfe_hits} / ${c.population} (${c.mfe_cov_pct}%)` : "&mdash;"}</td>
      <td><button type="button" class="btn btn-secondary probe-validate" data-seq="${escHtml(c.seq)}" ${running ? "disabled" : ""} title="이 Probe와 이 세트의 F/R로 실제 검증 과정(population·표준물질·dimer/hairpin)을 실행합니다. 결과는 이 카드에만 표시됩니다.">&#9654; 검증</button></td></tr>${valRow}`;
  }).join("");
  return `<p class="note">F/R 사이 정렬 열 ${r.col_window[0]}-${r.col_window[1]} 에서 ${r.n_pass_rules.toLocaleString()}개 창이 규칙을 통과했습니다
    (길이 ${r.len_range[0]}-${r.len_range[1]}nt · IUPAC 최대 ${r.max_iupac ?? 0}개 · GC ${r.gc_range[0]}-${r.gc_range[1]}% · Tm ${r.tm_range[0]}-${r.tm_range[1]}&deg;C = 더 높은 프라이머(${r.primer_tm}&deg;C)보다 ${r.tm_above ? `${r.tm_above[0]}-${r.tm_above[1]}` : "3-10"}&deg;C 높게 · 5′ G 제외 · 4개 이상 연속 같은 염기 제외).
    ${escHtml(r.generated)} 실행. ${pdScopeBadge(r.host_scope)} <span class="note">(집계 대상 ${r.population_n ?? "전체"}개)</span></p>
    <p class="note"><b>등록된 Probe</b> <span class="mono">${escHtml(reg.seq || "")}</span> · Tm ${reg.tm != null ? reg.tm : "&mdash;"}&deg;C · 정렬에서 완전 일치 ${reg.msa_cov0}${pct(reg.msa_cov0)}</p>
    ${rows ? `<div style="overflow:auto;"><table class="mini-reference-table"><tr><th>#</th><th>서열 (5′→3′)</th><th>출처</th><th>가닥</th><th>길이</th><th>IUPAC</th><th>Tm&deg;C <span class="note">(Primer3 Tm)</span></th><th>GC</th><th>정렬 열</th>
      <th>정렬 완전 일치(0 mm)</th><th>1 mm 이하</th><th>MFEprimer population (${r.probe_mismatch} mm)</th><th>실제 검증</th></tr>${rows}</table></div>
      <p class="note">Tm은 출처마다 자기 계산식 값이고(자체 규칙: nearest-neighbour, Primer3: Primer3 계산) 괄호는 같은 조건(250 nM, 50 mM Na⁺)의 Primer3 Tm이라 MFEprimer의 Tm과 몇 도 다를 수 있습니다. Primer3 후보는 Primer3의 internal oligo 설계(Tm·GC·poly-X·hairpin·penalty)로 뽑고 5′ G인 가닥은 뺐으며, IUPAC는 뒤에서 같은 방식으로 넣었습니다. 합성 전에 hairpin/dimer와 형광체 위치(5′ 형광체, 3′ 소광체)를 따로 확인하세요. 이 결과는 in silico 후보입니다.</p>`
      : "<p class='note'>규칙을 통과한 후보가 없습니다 (Tm/GC 조건이 이 영역에는 너무 엄격할 수 있습니다).</p>"}`;
}

function probeValidationHtml(v, setData, open) {
  const reg = setData.summary || {};
  const pctf = (x) => (x == null || x === "" ? "&mdash;" : `${Number(x).toFixed(1)}%`);
    const sm = v.summary || {};
    const delta = (a, b) => (a == null || b == null ? "" : ` <span class="note">(등록 Probe 세트 ${pctf(b)})</span>`);
    const std = v.standard || {};
    const dh = v.dimer_hairpin || {};
    const flagged = [...Object.values(dh.dimer || {}), ...Object.values(dh.hairpin || {})].reduce((n, x) => n + ((x && x.flagged) || 0), 0);
    const lin = (v.lineages || []).map((l) => `<tr><td>${escHtml(l.lineage)}</td><td>${l.validated} / ${l.total}</td><td>${l.total ? (100 * l.validated / l.total).toFixed(1) : "0.0"}%</td></tr>`).join("");
    return `<details class="probe-validation" ${open ? "open" : ""} style="margin-top:12px;border:1px solid var(--border);border-radius:var(--r-md);padding:8px 12px;">
      <summary><b class="mono">${escHtml(v.seq)}</b> <span class="note">· ${escHtml(v.generated)} 검증 결과</span> ${pdScopeBadge(v.host_scope)} <button type="button" class="btn btn-secondary probe-msa-btn" data-seq="${escHtml(v.seq)}" style="padding:0 8px;" title="이 Probe로 만든 전체 게놈 정렬(MSA)을 창으로 봅니다">&#129516; MSA 보기</button></summary>
      <div style="display:flex;gap:18px;flex-wrap:wrap;margin:8px 0;">
        <div><div class="note">최종 Validated</div><b>${sm.validated ?? "&mdash;"} / ${sm.total_genomes ?? "&mdash;"} (${pctf(sm.validated_pct)})</b>${delta(sm.validated_pct, reg.validated_pct)}</div>
        <div><div class="note">Probe 완전 일치</div><b>${pctf(sm.probe_coverage_mm0_pct)}</b>${delta(sm.probe_coverage_mm0_pct, reg.probe_coverage_mm0_pct)}</div>
        <div><div class="note">Primer 커버리지 (mm ${v.primer_mismatch_pass})</div><b>${pctf(sm["primer_coverage_mm" + v.primer_mismatch_pass + "_pct"])}</b></div>
        <div><div class="note">Probe dropout</div><b>${sm.probe_dropout ?? "&mdash;"}</b>${reg.probe_dropout != null ? ` <span class="note">(등록 ${reg.probe_dropout})</span>` : ""}</div>
        <div><div class="note">표준물질(KDCA)</div><b>${std.ok ? "✅ 매치" : (std.matched && std.matched.length ? "△ 일부" : "✗ 불일치")}</b></div>
        <div><div class="note">dimer/hairpin 경고</div><b>${flagged}</b></div>
      </div>
      <table class="mini-reference-table"><tr><th>${escHtml((PROFILE.lineage && PROFILE.lineage.label) || "Lineage")}</th><th>검증 통과 / 전체</th><th>%</th></tr>${lin}</table>
      <p class="note">실제 파이프라인과 같은 검색·판정으로 계산했습니다 (F/R mismatch&le;${v.primer_mismatch_pass} + 3′ 끝 일치, Probe mismatch&le;${v.probe_mismatch}). 이 결과는 이 카드에만 저장되며 개요·요약·보고서에는 들어가지 않습니다.</p>
    </details>`;
}

function renderProbeValidations(items, setData) {
  if (!items.length) return "";
  const blocks = items.map((v, vi) => probeValidationHtml(v, setData, vi === 0)).join("");
  return `<h4 style="margin:16px 0 0;">이전 검증 기록 <span class="note">(현재 후보 목록에 없는 Probe · 이 카드 전용)</span></h4>${blocks}`;
}

const PD_SCOPE_LABEL = { human: "사람 유래", animal: "동물 유래", unknown: "숙주 미상", all: "전체" };
function pdScopeBadge(scope) {
  if (!scope) return "";
  const stale = state.hostScope && state.hostScope !== scope;
  return `<span class="pill ${stale ? "pill-warn" : "pill-ok"}" title="이 결과를 계산할 때 쓴 숙주 기준입니다${stale ? ` — 지금 화면의 기준(${PD_SCOPE_LABEL[state.hostScope] || state.hostScope})과 다르니 다시 실행하세요` : ""}">숙주 기준: ${PD_SCOPE_LABEL[scope] || scope}</span>`;
}

function renderProbeBest(b) {
  if (!b || !b.generated) return "";
  const pctf = (v) => (v == null ? "&mdash;" : `${Number(v).toFixed(1)}%`);
  const rows = (b.ranking || []).map((r, i) => `<tr${i === 0 ? ' style="background:rgba(84,162,74,.14);font-weight:600;"' : ""}>
    <td>${i + 1}</td><td class="mono">${escHtml(r.seq)} <button type="button" class="btn btn-secondary probe-copy" data-seq="${escHtml(r.seq)}" style="padding:0 6px;" title="서열 복사">&#128203;</button></td>
    <td>${({ own: "자체", primer3: "Primer3", "own+primer3": "자체+P3" })[r.source] || ""}</td><td>${r.n_iupac}</td><td>${r.tm ?? "&mdash;"}</td><td>${r.validated ?? "&mdash;"} (${pctf(r.validated_pct)})</td><td>${pctf(r.probe_pct)}</td><td>${r.dropout ?? "&mdash;"}</td>
    <td>${r.std_ok ? "✅" : "✗"}</td></tr>`).join("");
  const top = (b.ranking || [])[0];
  return `<h4 style="margin:16px 0 0;">&#127942; 최적 Probe 탐색 결과 <span class="note">${escHtml(b.generated)} · 후보 ${b.n_tested}개를 모두 실제 검증</span> ${pdScopeBadge(b.host_scope)}</h4>
    ${top ? `<p style="margin:6px 0;">가장 좋은 Probe: <b class="mono">${escHtml(top.seq)}</b> — Validated ${top.validated} (${pctf(top.validated_pct)}), Probe 완전 일치 ${pctf(top.probe_pct)}, dropout ${top.dropout}</p>` : ""}
    <div style="overflow:auto;"><table class="mini-reference-table"><tr><th>순위</th><th>서열 (5′→3′)</th><th>출처</th><th>IUPAC</th><th>Tm&deg;C</th><th>최종 Validated</th><th>Probe 완전 일치</th><th>dropout</th><th>표준물질</th></tr>${rows}</table></div>
    <p class="note">순위 기준: 최종 Validated 수 &rarr; 표준물질 매치 &rarr; Probe 완전 일치 &rarr; dropout이 적은 순 &rarr; IUPAC 적은 순. 후보마다 자세한 검증 결과는 위 후보 표의 각 행 아래에 있습니다. 최종 Validated는 F/R이 정하는 상한이 있어서, Probe만으로는 차이가 작을 수 있습니다. 이 결과도 이 카드에만 저장됩니다.</p>`;
}

// The alignment of one designed probe (F/R of the set + that probe), in a modal -- the same view as the set's MSA card.
async function openProbeMsaModal(setData, seq, urlOverride, titleLabel) {
  let d;
  try {
    d = await fetchJSON(urlOverride || `/api/sets/${encodeURIComponent(setData.set_id)}/probe-design/msa?host=${encodeURIComponent(state.hostScope || "")}&seq=${encodeURIComponent(seq)}`);
  } catch (e) { uiAlert("MSA를 불러오지 못했습니다. 이 Probe의 검증 결과가 있는지 확인하세요."); return; }
  const pseudo = { ...d, set_id: `${setData.set_id} + ${seq}` };
  const saved = { s: state.msaSearch, m: state.msaFilterMode, c: state.msaColFilters, a: state._msaAttrOf, b: state._msaBase };
  state.msaSearch = ""; state.msaFilterMode = "all"; state.msaColFilters = {};
  const { overlay } = openDetailModal({
    title: titleLabel || `전체 게놈 정렬 (MSA) &mdash; 설계 Probe <span class="mono">${escHtml(seq)}</span>`,
    subtitle: titleLabel ? "varVAMP가 설계한 F/Probe/R · 숙주 기준을 따릅니다 · 컬럼은 ⌥/Alt+드래그로 선택" : `F/R은 ${escHtml(setData.display_name || setData.set_id)} 그대로, Probe만 설계한 서열입니다 · 숙주 기준을 따릅니다 · 컬럼은 ⌥/Alt+드래그로 선택`,
    html: '<div class="card msa-card" id="probeMsaCard" style="margin:0;"></div>',
    compact: true,
  });
  const box = overlay.querySelector(".modal-box");
  box.style.width = "min(97vw, 1500px)"; box.style.maxWidth = "97vw";
  const card = overlay.querySelector("#probeMsaCard");
  const HEADER = `<div class="msa-card-header"><span class="note">${d.msa && d.msa.rows ? d.msa.rows.length : 0}개 게놈</span>
    <div style="display:flex;gap:6px;flex-wrap:wrap;">
      <button type="button" class="btn btn-secondary pm-snap">&#128247; 이미지 저장</button>
      <button type="button" class="btn btn-secondary pm-copyimg">&#128203; 이미지 복사</button>
      <button type="button" class="btn btn-secondary pm-copytxt">&#128203; 텍스트 복사</button>
      <button type="button" class="btn btn-secondary pm-fs">&#9974; 전체화면</button></div></div>`;
  const draw = () => {
    const focusSearch = document.activeElement && document.activeElement.id === "msaSearch" && card.contains(document.activeElement);
    const caret = focusSearch ? document.activeElement.selectionStart : null;
    card.innerHTML = HEADER + renderMSA(pseudo);
    const pi = card.querySelector("#msaDiversityRow");              // the pi track belongs to the registered set only
    if (pi) pi.remove();
    const ruler = card.querySelector(".msa-ruler");
    if (ruler) ruler.style.top = "0";
    const comb = card.querySelector(".msa-combined-row"), cons = card.querySelector(".msa-consensus");      // the sticky rows move up with the pi track gone
    if (comb) comb.style.top = "17px";
    if (cons) cons.style.top = "35px";
    const q = (sel) => card.querySelector(sel);
    q("#msaSearch")?.addEventListener("input", (e) => { state.msaSearch = e.target.value; draw(); });
    card.querySelectorAll(".mt-funnel[data-msa-col]").forEach((b) =>
      b.addEventListener("click", (e) => { e.stopPropagation(); msaOpenColMenu(b, b.dataset.msaCol, draw); }));
    q("#msaClearColFilters")?.addEventListener("click", () => { state.msaColFilters = {}; draw(); });
    q("#msaMetaCols")?.addEventListener("change", (e) => { state.msaMetaCols = e.target.checked; draw(); });
    q("#msaFilterMode")?.addEventListener("change", (e) => { state.msaFilterMode = e.target.value; draw(); });
    q(".pm-snap")?.addEventListener("click", () => msaSnapshot(card, pseudo, "save"));
    q(".pm-copyimg")?.addEventListener("click", () => msaSnapshot(card, pseudo, "copy"));
    q(".pm-copytxt")?.addEventListener("click", () => msaCopyText(card, pseudo));
    q(".pm-fs")?.addEventListener("click", () => toggleFullscreen(card));
    const copyBtn = q("#msaColSelectCopyBtn"), clearBtn = q("#msaColSelectClearBtn");
    if (copyBtn) copyBtn.addEventListener("click", () => msaCopyColSelection(card));
    if (clearBtn) clearBtn.addEventListener("click", () => msaClearColSelection(card));
    if (focusSearch) { const el2 = q("#msaSearch"); el2.focus(); el2.setSelectionRange(caret, caret); }
  };
  draw();
  wireColumnSelection(card);
  const mo = new MutationObserver(() => {                       // closed: give the set page its own filter state back
    if (document.body.contains(overlay)) return;
    state.msaSearch = saved.s; state.msaFilterMode = saved.m; state.msaColFilters = saved.c; state._msaAttrOf = saved.a; state._msaBase = saved.b;
    mo.disconnect();
  });
  mo.observe(document.body, { childList: true, subtree: true });
}

async function initProbeDesign(detail, setData) {
  const body = detail.querySelector("#probeDesignBody"), btn = detail.querySelector("#probeDesignRunBtn");
  if (!body || !btn) return;
  const id = encodeURIComponent(setData.set_id);
  const valBox = detail.querySelector("#probeDesignValidations");
  let design = null, vals = [], running = null;       // the last design, the saved validations, the validation in progress
  const rerender = () => {
    if (!document.body.contains(body)) return;
    body.innerHTML = renderProbeDesign(design, { vals: new Map(vals.map((v) => [v.seq, v])), running, setData });
    const inTable = new Set(((design && design.candidates) || []).map((c) => c.seq));
    valBox.innerHTML = renderProbeValidations(vals.filter((v) => !inTable.has(v.seq)), setData);
  };
  const show = (r) => { design = r; rerender(); };
  const loadValidations = async () => {
    try { vals = await fetchJSON(`/api/sets/${id}/probe-design/validations`); } catch (e) { /* none yet */ }
    rerender();
  };
  try { design = await fetchJSON(`/api/sets/${id}/probe-design`); } catch (e) { body.innerHTML = "<p class='note'>결과를 불러오지 못했습니다.</p>"; }
  await loadValidations();
  const bestBox = detail.querySelector("#probeDesignBest");
  const loadBest = async (status) => {
    let b = null;
    try { b = await fetchJSON(`/api/sets/${id}/probe-design/best`); } catch (e) { /* none yet */ }
    if (document.body.contains(bestBox)) bestBox.innerHTML = (status ? `<p class="note">${status}</p>` : "") + renderProbeBest(b);
  };
  loadBest();
  detail.querySelector("#probeBestBtn").addEventListener("click", async (e) => {
    const bb = e.currentTarget;
    const n = Number(detail.querySelector("#probeDesignTop").value) || 20;
    if (!(await uiConfirm(`후보 ${n}개를 모두 실제 검증 과정에 태웁니다. 약 ${Math.ceil(n * 7 / 60)}분 걸리며 그동안 서버가 계산에 쓰입니다. 계속할까요?`, { title: "최적 Probe 자동 탐색", okLabel: "시작" }))) return;
    bb.disabled = true;
    try {
      const { job_id, error } = await fetchJSON(`/api/sets/${id}/probe-design/best`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ top: n, max_iupac: Number(detail.querySelector("#probeDesignIupac").value), method: detail.querySelector("#probeDesignMethod").value, host: state.hostScope,
          tm_above_min: Number(detail.querySelector("#probeDesignTmMin").value), tm_above_max: Number(detail.querySelector("#probeDesignTmMax").value),
          probe_mismatch: Number(detail.querySelector("#probeDesignMm").value) || 0 }) });
      if (error) throw new Error(error);
      const jt = jobToast("Probe 최적 탐색");
      for (;;) {
        const j = await fetchJSON(`/api/jobs/${job_id}`);
        if (!document.body.contains(bestBox)) { jt.done(j.status !== "error"); return; }
        if (j.status === "running") { jt.progress(j); bestBox.querySelector(".probe-best-status")?.remove();
          bestBox.insertAdjacentHTML("afterbegin", `<p class="note probe-best-status">⏳ ${escHtml(j.progress || "")} (${j.elapsed_sec ?? 0}초)</p>`); await new Promise((r) => setTimeout(r, 2000)); continue; }
        if (j.status === "error") uiAlert("최적 탐색에 실패했습니다: " + (j.error || ""));
        jt.done(j.status !== "error");
        break;
      }
      await Promise.all([loadBest(), loadValidations()]);
    } catch (err) { uiAlert("최적 탐색을 시작하지 못했습니다: " + err.message); }
    bb.disabled = false;
  });
  detail.querySelector("#probeDesignCard").addEventListener("click", async (e) => {
    const mb = e.target.closest(".probe-msa-btn");
    if (mb) { e.preventDefault(); openProbeMsaModal(setData, mb.dataset.seq); return; }
    const vb = e.target.closest(".probe-validate");
    if (vb) {
      if (running) { uiAlert("이미 다른 Probe를 검증 중입니다. 끝난 뒤 다시 눌러 주세요."); return; }
      running = { seq: vb.dataset.seq, text: "시작..." }; rerender();
      try {
        const { job_id, error } = await fetchJSON(`/api/sets/${id}/probe-design/validate`, {
          method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ seq: vb.dataset.seq, host: state.hostScope }) });
        if (error) throw new Error(error);
        const jt = jobToast("Probe 검증");
        for (;;) {
          const j = await fetchJSON(`/api/jobs/${job_id}`);
          if (!document.body.contains(valBox)) { jt.done(j.status !== "error"); return; }
          if (j.status === "running") { jt.progress(j); running.text = j.progress || ""; rerender(); await new Promise((r) => setTimeout(r, 1500)); continue; }
          if (j.status === "error") uiAlert("검증에 실패했습니다: " + (j.error || ""));
          jt.done(j.status !== "error");
          break;
        }
      } catch (err) { uiAlert("검증을 시작하지 못했습니다: " + err.message); }
      running = null;
      await loadValidations();
      return;
    }
    const c = e.target.closest(".probe-copy");
    if (c) copyTextToClipboard(c.dataset.seq).then(() => { c.textContent = "✓"; setTimeout(() => { c.innerHTML = "&#128203;"; }, 1000); }).catch(() => uiAlert("복사하지 못했습니다."));
  });
  btn.addEventListener("click", async () => {
    btn.disabled = true; const orig = btn.innerHTML; btn.textContent = "설계·검증 중...";
    toast("Probe 설계·검증 중...", { type: "run", id: "probe-design" });
    try {
      const r = await fetchJSON(`/api/sets/${id}/probe-design`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ top: Number(detail.querySelector("#probeDesignTop").value) || 10, max_iupac: Number(detail.querySelector("#probeDesignIupac").value), method: detail.querySelector("#probeDesignMethod").value, host: state.hostScope, tm_above_min: Number(detail.querySelector("#probeDesignTmMin").value), tm_above_max: Number(detail.querySelector("#probeDesignTmMax").value), probe_mismatch: Number(detail.querySelector("#probeDesignMm").value) || 0 }),
      });
      if (r.error) { uiAlert(r.error); toast("Probe 설계 실패", { type: "err", id: "probe-design" }); } else { show(r); toast("Probe 설계 완료", { type: "ok", id: "probe-design" }); }
    } catch (e) { uiAlert("Probe 설계에 실패했습니다: " + e.message); toast("Probe 설계 실패", { type: "err", id: "probe-design" }); }
    btn.disabled = false; btn.innerHTML = orig;
  });
}

// navigator.clipboard only exists on HTTPS/localhost; the dashboard is usually opened over plain http on the LAN
// (http://10.x.x.x:8765), so text copy falls back to a hidden textarea + execCommand("copy").
function copyTextToClipboard(text) {
  const legacy = () => new Promise((resolve, reject) => {
    const ta = document.createElement("textarea");
    ta.value = text;
    ta.style.cssText = "position:fixed;left:-9999px;top:0;opacity:0;";
    document.body.appendChild(ta);
    ta.select();
    let ok = false;
    try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
    ta.remove();
    ok ? resolve() : reject(new Error("copy failed"));
  });
  // the async API can refuse (page not focused, permission): fall back to the old way before giving up
  if (navigator.clipboard && window.isSecureContext) return navigator.clipboard.writeText(text).catch(legacy);
  return legacy();
}

// the column window to export: the dragged selection if there is one, otherwise the whole alignment
function msaExportWindow(msaCard, setData) {
  const n = setData.msa.columns;
  const readout = msaCard.querySelector("#msaColSelectReadout");
  if (readout && !readout.hidden && readout.dataset.start !== undefined) {
    const lo = Math.max(0, Number(readout.dataset.start)), hi = Math.min(n - 1, Number(readout.dataset.end));
    if (!isNaN(lo) && !isNaN(hi) && hi >= lo) return { lo, hi, selected: true };
  }
  return { lo: 0, hi: n - 1, selected: false };
}

// every row that is currently listed (search/filters applied, NOT limited to the 800 drawn), as {label, seq, kind}
function msaExportRows(msaCard, setData) {
  const msa = setData.msa;
  const textOf = (sel) => { const el = msaCard.querySelector(sel); return el ? el.textContent : ""; };
  const out = [
    { kind: "ruler", label: "ref position", seq: textOf(".msa-ruler .msa-seq") },
    { kind: "combined", label: "PRIMER/PROBE", seq: textOf(".msa-combined-row .msa-seq") },
    { kind: "consensus", label: "CONSENSUS", seq: msa.consensus },
  ];
  for (const r of ((state._msaRowsBySet && state._msaRowsBySet.get(setData)) || msa.rows)) {
    out.push({
      kind: "genome",
      label: `${r.genome_id} | ${r.lineage || "NA"} | ${r.country || "?"}  F${mmText(r.f_mismatch)} P${mmText(r.p_mismatch)} R${mmText(r.r_mismatch)}`,
      seq: r.dots, tsvLabel: `${r.genome_id} | ${r.lineage || "NA"} | ${r.country || "?"}`,
    });
  }
  return out;
}

function msaCopyText(msaCard, setData) {
  const msa = setData.msa;
  const { lo, hi, selected } = msaExportWindow(msaCard, setData);
  const undot = (d) => d.split("").map((c, i) => (c === "." ? msa.consensus[i] : c)).join("");
  const lines = msaExportRows(msaCard, setData).map((r) => {
    const seq = r.kind === "genome" ? undot(r.seq) : r.seq;
    return `${r.tsvLabel || r.label}\t${seq.slice(lo, hi + 1)}`;
  });
  copyTextToClipboard(lines.join("\n") + "\n")
    .then(() => uiAlert(`정렬 텍스트를 복사했습니다 (${lines.length - 3}개 게놈, 컬럼 ${lo + 1}-${hi + 1}${selected ? ", 선택 구간" : ", 전체"}).`))
    .catch(() => uiAlert("클립보드 복사에 실패했습니다. 브라우저가 복사를 막았을 수 있습니다."));
}

// Draws the alignment on a canvas (no external library): the three header rows + every listed genome.
// "." (same as consensus) stays faint, a letter (mismatch) is red, "-" is grey; F/Probe/R windows are shaded like on screen.
function msaRenderCanvas(msaCard, setData) {
  const msa = setData.msa;
  const { lo, hi } = msaExportWindow(msaCard, setData);
  const rows = msaExportRows(msaCard, setData);
  const cs = getComputedStyle(document.documentElement);
  const css = (v, d) => (cs.getPropertyValue(v).trim() || d);
  const bg = css("--bg", "#ffffff"), fg = css("--fg", "#1a1d23"), muted = css("--muted", "#6b7280"), border = css("--border", "#e5e7eb");
  const FS = 12, LH = 16, PAD = 8, LABEL_W = 330;
  const font = `${FS}px ui-monospace, SFMono-Regular, Menlo, Consolas, monospace`;
  const probe = document.createElement("canvas").getContext("2d");
  probe.font = font;
  const cw = probe.measureText("M").width;
  const ncol = hi - lo + 1;
  const W = Math.ceil(PAD * 2 + LABEL_W + ncol * cw), H = PAD * 2 + rows.length * LH;
  // browsers cap canvas size (~32k px per side, ~270M px total): shrink the drawing if the alignment is huge
  const scale = Math.min(1, 30000 / W, 30000 / H, Math.sqrt(120e6 / (W * H)));
  const canvas = document.createElement("canvas");
  canvas.width = Math.max(1, Math.floor(W * scale));
  canvas.height = Math.max(1, Math.floor(H * scale));
  const ctx = canvas.getContext("2d");
  ctx.scale(scale, scale);
  ctx.fillStyle = bg; ctx.fillRect(0, 0, W, H);
  ctx.font = font; ctx.textBaseline = "middle";
  const wins = [[msa.f_col, "rgba(76,120,168,.32)"], [msa.p_col, "rgba(242,166,90,.42)"], [msa.r_col, "rgba(84,162,74,.32)"]];
  const x0 = PAD + LABEL_W;
  rows.forEach((r, ri) => {
    const y = PAD + ri * LH;
    if (r.kind !== "ruler") for (const [w, color] of wins) {
      if (!w) continue;
      const a = Math.max(w[0], lo), b = Math.min(w[1], hi);
      if (b >= a) { ctx.fillStyle = color; ctx.fillRect(x0 + (a - lo) * cw, y, (b - a + 1) * cw, LH); }
    }
    ctx.fillStyle = r.kind === "genome" ? fg : muted;
    ctx.font = (r.kind === "genome" ? "" : "bold ") + font;
    let label = r.label;
    while (label.length > 4 && ctx.measureText(label).width > LABEL_W - 8) label = label.slice(0, -2);
    if (label !== r.label) label += "…";
    ctx.fillText(label, PAD, y + LH / 2);
    ctx.font = font;
    const seq = r.seq.slice(lo, hi + 1);
    for (let i = 0; i < seq.length; i++) {
      const c = seq[i];
      if (c === " ") continue;
      if (r.kind === "genome") ctx.fillStyle = c === "." ? border : c === "-" ? muted : "#E45756";
      else ctx.fillStyle = fg;
      ctx.fillText(c, x0 + i * cw, y + LH / 2);
    }
  });
  ctx.strokeStyle = border; ctx.beginPath(); ctx.moveTo(PAD + LABEL_W - 4, PAD); ctx.lineTo(PAD + LABEL_W - 4, H - PAD); ctx.stroke();
  return canvas;
}

function msaSnapshot(msaCard, setData, mode) {
  if (!setData.msa || !setData.msa.columns) return;
  const canvas = msaRenderCanvas(msaCard, setData);
  const name = `${setData.set_id}_alignment.png`;
  const download = (blob) => {
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url; a.download = name;
    document.body.appendChild(a); a.click(); a.remove();
    setTimeout(() => URL.revokeObjectURL(url), 2000);
  };
  canvas.toBlob((blob) => {
    if (!blob) { uiAlert("이미지를 만들지 못했습니다. 정렬이 너무 큽니다 — 컬럼을 드래그로 선택해 구간만 저장해 보세요."); return; }
    if (mode === "copy") {
      if (navigator.clipboard && window.ClipboardItem && window.isSecureContext) {
        navigator.clipboard.write([new ClipboardItem({ "image/png": blob })])
          .then(() => uiAlert("정렬 이미지를 클립보드에 복사했습니다."))
          .catch(() => { download(blob); uiAlert("클립보드 복사가 막혀 PNG 파일로 저장했습니다."); });
      } else {
        download(blob);
        uiAlert("이 주소(http)에서는 브라우저가 이미지 복사를 허용하지 않아 PNG 파일로 저장했습니다. (localhost 또는 HTTPS에서는 바로 복사됩니다.)");
      }
    } else download(blob);
  }, "image/png");
}

function msaCopyColSelection(msaCard) {
  const readout = msaCard.querySelector("#msaColSelectReadout");
  if (!readout || readout.hidden) return;
  const lo = Number(readout.dataset.start), hi = Number(readout.dataset.end);
  const lines = [];
  msaCard.querySelectorAll(".msa-row").forEach((row) => {
    const label = row.dataset.label || row.querySelector(".msa-label")?.textContent?.trim().replace(/\s+/g, " ") || "";
    const seqEl = row.querySelector(".msa-seq");
    if (!seqEl) return;
    const text = seqEl.textContent;
    lines.push(`${label}\t${text.slice(lo, hi + 1)}`);
  });
  const text = lines.join("\n");
  copyTextToClipboard(text).then(() => {
    const btn = msaCard.querySelector("#msaColSelectCopyBtn");
    if (btn) {
      const orig = btn.textContent;
      btn.textContent = "복사됨!";
      setTimeout(() => { btn.textContent = orig; }, 1200);
    }
  }).catch(() => uiAlert("클립보드 복사에 실패했습니다."));
}

function wireColumnSelection(msaCard) {
  msaCard.addEventListener("mousedown", (e) => {
    if (!e.altKey) return; // VS Code-style: Alt+drag enters column/box-select mode
    const outer = e.target.closest(".msa-outer");
    if (!outer || e.target.closest(".msa-label")) return;
    e.preventDefault();
    const col = msaColAt(msaCard, e.clientX);
    msaColSelectDrag = { msaCard, startCol: col };
    msaUpdateColOverlay(msaCard, col, col);
  });
  const copyBtn = msaCard.querySelector("#msaColSelectCopyBtn");
  if (copyBtn) copyBtn.addEventListener("click", () => msaCopyColSelection(msaCard));
  const clearBtn = msaCard.querySelector("#msaColSelectClearBtn");
  if (clearBtn) clearBtn.addEventListener("click", () => msaClearColSelection(msaCard));
}

if (!window.__msaColSelectGlobalWired) {
  window.__msaColSelectGlobalWired = true;
  document.addEventListener("mousemove", (e) => {
    if (!msaColSelectDrag) return;
    const { msaCard, startCol } = msaColSelectDrag;
    if (!document.body.contains(msaCard)) { msaColSelectDrag = null; return; }
    const col = msaColAt(msaCard, e.clientX);
    msaUpdateColOverlay(msaCard, startCol, col);
  });
  document.addEventListener("mouseup", (e) => {
    if (!msaColSelectDrag) return;
    const { msaCard, startCol } = msaColSelectDrag;
    if (document.body.contains(msaCard)) {
      const col = msaColAt(msaCard, e.clientX);
      msaShowColReadout(msaCard, startCol, col);
    }
    msaColSelectDrag = null;
  });
}

// Binding detail (primer/probe alignment on this genome's amplicon) in a modal window -- opened
// from a row of the per-genome tables and from a genome row of the MSA.
function openGenomeModal(g, setData, labels, combo) {
  const status = comboRowStatus(g);
  const meta = [g.lineage || "NA", g.country || "", combo || ""].filter(Boolean).map(escHtml).join(" · ");
  return openDetailModal({
    title: `<span class="mono">${escHtml(g.genome_id)}</span> <span class="combo-status combo-status-${status.level}">${status.icon} ${status.label}</span>`,
    subtitle: meta,
    html: highlightAmplicon(g, setData, labels),
  });
}

function renderBinding(binding) {
  if (!binding || !binding.length) return "<div class='note'>참조유전체 결합 정보 없음</div>";
  return `<div class="binding-track">${binding.map((h) => {
    if (!h.found) {
      return `<div class="binding-row"><span class="role ${h.role}">${h.role}</span><span style="color:var(--bad);font-weight:600;">NO MATCH (참조유전체에서 결합 못 찾음)</span></div>`;
    }
    return `<div class="binding-row"><span class="role ${h.role}">${h.role}</span>
      <span>${h.chrom}:${h.start}-${h.end} (${h.strand}strand)</span>
      <span class="mm ${mmClass(h.mismatches)}" style="padding:1px 8px;border-radius:5px;">mismatch ${h.mismatches}</span></div>`;
  }).join("")}</div>`;
}

const DIMER_PAIR_LABEL = {
  "F-F": "Forward self-dimer", "P-P": "Probe self-dimer", "R-R": "Reverse self-dimer",
  "F-P": "Forward &harr; Probe", "F-R": "Forward &harr; Reverse", "P-R": "Probe &harr; Reverse",
};
const DIMER_ROLE_LABEL = { F: "Forward", P: "Probe", R: "Reverse" };

function renderDimerHairpin(dh) {
  if (!dh) return "<div class='note'>계산되지 않음</div>";
  const threshold = dh.dg_flag_threshold ?? -9;
  const dimerEntries = Object.entries(dh.dimer || {});
  const hairpinEntries = Object.entries(dh.hairpin || {});
  if (!dimerEntries.length && !hairpinEntries.length) {
    return `<div class="note">감지된 dimer/hairpin 없음 (mfeprimer score cutoff 이상)</div>
      <div class="note" style="margin-top:4px;">&Delta;G &le; ${threshold} kcal/mol이면 &ldquo;주의&rdquo;로 표시하는 참고 기준선입니다.</div>`;
  }
  const row = (label, dg, extraHtml, align) => {
    const bad = dg <= threshold;
    return `<div class="dimer-row ${bad ? "bad" : ""}" title="${escHtml(align || "")}">
      <span class="dimer-label">${label}</span>
      <span class="dimer-dg mono">&#916;G ${dg.toFixed(2)}</span>
      ${extraHtml || ""}
      ${bad ? `<span class="pill pill-bad" style="font-size:9px;">주의</span>` : `<span class="pill pill-ok" style="font-size:9px;">OK</span>`}
    </div>`;
  };
  const dimerHtml = dimerEntries
    .sort((a, b) => a[1].worst_dg - b[1].worst_dg)
    .map(([pair, v]) => row(DIMER_PAIR_LABEL[pair] || pair, v.worst_dg, null, v.align))
    .join("");
  const hairpinHtml = hairpinEntries
    .sort((a, b) => a[1].worst_dg - b[1].worst_dg)
    .map(([role, v]) => row(
      `${DIMER_ROLE_LABEL[role] || role} hairpin`, v.worst_dg,
      `<span class="dimer-tm mono">Tm ${v.worst_tm.toFixed(1)}&deg;C</span>`, v.align,
    ))
    .join("");
  return `<div class="dimer-hairpin-list">${dimerHtml}${hairpinHtml}</div>
    <div class="note" style="margin-top:6px;">&Delta;G &le; ${threshold} kcal/mol이면 &ldquo;주의&rdquo; &mdash; 결합/헤어핀 형성으로 증폭 효율이 떨어질 수 있다는 참고 기준선입니다(IDT OligoAnalyzer 등에서 흔히 쓰는 경험적 값). 각 항목에 마우스를 올리면 실제 결합 형태를 볼 수 있습니다. IUPAC 축퇴염기가 있는 프라이머는 모든 조합 중 가장 강한(가장 음수인) 경우를 보여줍니다.</div>`;
}

function markRanges(seq, ranges) {
  // ranges: [{start, end, cls}], non-overlapping, sorted by start. Wraps the given
  // [start,end) slices in <mark class="cls">, leaving everything else untouched.
  const sorted = ranges.filter((r) => r && r.start != null).sort((a, b) => a.start - b.start);
  let out = "";
  let i = 0;
  for (const r of sorted) {
    if (r.start > i) out += seq.slice(i, r.start);
    out += `<mark class="${r.cls}">${seq.slice(r.start, r.end)}</mark>`;
    i = r.end;
  }
  out += seq.slice(i);
  return out;
}

const IUPAC_JS = {
  A: "A", C: "C", G: "G", T: "T",
  R: "AG", Y: "CT", S: "GC", W: "AT", K: "GT", M: "AC",
  B: "CGT", D: "AGT", H: "ACT", V: "ACG", N: "ACGT",
};
function iupacMatches(oligoChar, consensusChar) {
  if (oligoChar === "-" || consensusChar === "-") return true;
  const set = IUPAC_JS[oligoChar.toUpperCase()] || oligoChar.toUpperCase();
  return set.includes(consensusChar.toUpperCase());
}

function buildCombinedSeqRow(msa) {
  // One combined line spanning FORWARD/PROBE/REVERSE, showing each oligo's own literal
  // sequence letters (not an arrow decoration) so it can be read directly against the
  // CONSENSUS row right below it. A position where the oligo's IUPAC code doesn't cover
  // the population consensus base is flagged as a mismatch (extra "mm" class) -- this is
  // the primer/probe's own mismatch against the population, distinct from the per-genome
  // dot-diff rows below. Outside the three oligo spans it falls back to that oligo's own
  // aligned "-" gap text (from the shared mafft run), never blank space, so the line has
  // the same length/columns as every other alignment row.
  const n = msa.columns;
  const chars = new Array(n).fill("-");
  const cls = new Array(n).fill(null);
  const specs = [
    { key: "forward", c: "f" },
    { key: "probe", c: "p" },
    { key: "reverse", c: "r" },
  ];
  for (const { key, c } of specs) {
    const rr = msa.ref_rows && msa.ref_rows[key];
    if (!rr || !rr.col || !rr.aligned) continue;
    const [a, b] = rr.col;
    for (let i = a; i < b && i < n; i++) {
      const ch = rr.aligned[i];
      if (ch === "-") continue;
      chars[i] = ch;
      cls[i] = iupacMatches(ch, msa.consensus[i]) ? c : `${c} mm`;
    }
  }
  return { text: chars.join(""), cls };
}

function markByClsArray(text, clsArr) {
  // Wraps contiguous runs of the same non-null cls[i] in <mark class="cls">, leaving
  // unclassed characters (gap filler) untouched -- like markRanges() but driven by a
  // per-character class array instead of a list of [start,end) ranges.
  let out = "";
  let i = 0;
  while (i < text.length) {
    const c = clsArr[i];
    let j = i + 1;
    while (j < text.length && clsArr[j] === c) j++;
    const chunk = text.slice(i, j);
    out += c ? `<mark class="${c}">${chunk}</mark>` : chunk;
    i = j;
  }
  return out;
}

function strandLabel(strand) {
  if (strand === "+") return "정방향";
  if (strand === "-") return "역방향";
  return "";
}

// Value-filter menu for one MSA column (continent / country / lineage / verdict): checkboxes with counts, search, all/none.
function msaOpenColMenu(anchor, key, refresh) {
  document.querySelectorAll(".msa-col-menu").forEach((m) => m.remove());
  const attrOf = state._msaAttrOf, base = state._msaBase || [];
  const others = Object.entries(state.msaColFilters || {}).filter(([k]) => k !== key);
  const counts = new Map();
  for (const r of base) {
    const a = attrOf(r);
    if (!others.every(([k, set]) => set.has(a[k]))) continue;         // counts follow the other columns' filters
    counts.set(a[key], (counts.get(a[key]) || 0) + 1);
  }
  const all = [...counts.keys()].sort((a, b) => (a === "" ? 1 : b === "" ? -1 : a.localeCompare(b, undefined, { numeric: true })));
  const selected = new Set((state.msaColFilters || {})[key] || all);
  const menu = el("div", "mt-menu msa-col-menu");
  const rect = anchor.getBoundingClientRect();
  menu.style.left = `${Math.max(8, Math.min(rect.left, window.innerWidth - 270))}px`;
  menu.style.top = `${rect.bottom + 4}px`;
  const title = { region: "대륙", country: "나라", lineage: (PROFILE.lineage && PROFILE.lineage.label) || "Lineage", bind: "결합 여부", host: "숙주" }[key];
  menu.innerHTML = `<div class="mt-menu-head">${escHtml(title)} 필터</div>
    <input type="text" class="mt-menu-search" placeholder="값 검색..." spellcheck="false">
    <div class="mt-menu-actions"><a href="#" data-act="all">모두 선택</a><a href="#" data-act="none">모두 해제</a></div>
    <div class="mt-menu-list"></div>`;
  modalHost().appendChild(menu);
  const listEl = menu.querySelector(".mt-menu-list"), search = menu.querySelector(".mt-menu-search");
  const visible = () => { const q = search.value.trim().toLowerCase(); return all.filter((v) => !q || (v === "" ? "(비어 있음)" : v).toLowerCase().includes(q)); };
  const draw = () => {
    listEl.innerHTML = visible().slice(0, 300).map((v) => `<label class="mt-menu-item"><input type="checkbox" data-v="${escHtml(v)}" ${selected.has(v) ? "checked" : ""}>
      <span class="mt-menu-val">${v === "" ? "<i>(비어 있음)</i>" : escHtml(v)}</span><span class="mt-menu-n">${counts.get(v)}</span></label>`).join("") || "<div class='note' style='padding:6px 8px;'>일치하는 값 없음</div>";
  };
  const close = () => { menu.remove(); document.removeEventListener("mousedown", outside, true); document.removeEventListener("keydown", onKey, true); };
  const outside = (ev) => { if (!menu.contains(ev.target)) close(); };
  const onKey = (ev) => { if (ev.key === "Escape") { ev.stopPropagation(); close(); } };
  const apply = () => {
    state.msaColFilters = state.msaColFilters || {};
    if (selected.size === all.length) delete state.msaColFilters[key]; else state.msaColFilters[key] = new Set(selected);
    close(); refresh();                                              // re-renders the card (the menu is re-opened by the click)
    setTimeout(() => { const b = document.querySelector(`.mt-funnel[data-msa-col="${key}"]`); if (b) msaOpenColMenu(b, key, refresh); }, 0);
  };
  listEl.addEventListener("change", (ev) => {
    const cb = ev.target.closest("input[data-v]");
    if (!cb) return;
    if (cb.checked) selected.add(cb.dataset.v); else selected.delete(cb.dataset.v);
    apply();
  });
  search.addEventListener("input", draw);
  menu.querySelectorAll(".mt-menu-actions a").forEach((a) => a.addEventListener("click", (ev) => {
    ev.preventDefault();
    visible().forEach((v) => (a.dataset.act === "all" ? selected.add(v) : selected.delete(v)));
    apply();
  }));
  draw(); search.focus();
  document.addEventListener("mousedown", outside, true);
  document.addEventListener("keydown", onKey, true);
}

function renderMSA(setData) {
  const msa = setData.msa;
  if (!msa || !msa.columns) {
    return "<p class='note'>이 세트는 population에서 매치되는 게놈이 없어 정렬을 만들 수 없습니다.</p>";
  }
  const ranges = [
    msa.f_col && { start: msa.f_col[0], end: msa.f_col[1], cls: "f" },
    msa.p_col && { start: msa.p_col[0], end: msa.p_col[1], cls: "p" },
    msa.r_col && { start: msa.r_col[0], end: msa.r_col[1], cls: "r" },
  ].filter(Boolean);

  const search = (state.msaSearch || "").toLowerCase();
  let rows = msa.rows;
  const mode = state.msaFilterMode || "all";
  // One verdict per genome, the same one the genome table / pop-ups use (mismatches + binding quality): 정상, 주의 or 미결합.
  const meta = new Map((setData.genomes || []).map((g) => [g.genome_id, g]));
  const statusOf = (r) => comboRowStatus(meta.get(r.genome_id) || r);
  const bindState = (r) => { const lv = statusOf(r).level; return lv === "ok" ? "ok" : lv === "warn" ? "warn" : "bad"; };
  // the information columns are shown by default; the toolbar switch hides them to give the alignment more width
  const metaCols = state.msaMetaCols === undefined ? true : state.msaMetaCols;
  const funnel = (key) => `<button type="button" class="mt-funnel${(state.msaColFilters || {})[key] ? " on" : ""}" data-msa-col="${key}" title="이 열의 값으로 거르기">&#9662;</button>`;
  const counts = { ok: 0, warn: 0, bad: 0 };
  for (const r of msa.rows) counts[bindState(r)] += 1;
  if (mode === "bound") rows = rows.filter((r) => bindState(r) !== "bad");
  else if (mode === "ok" || mode === "warn" || mode === "bad") rows = rows.filter((r) => bindState(r) === mode);
  if (search) {
    rows = rows.filter((r) => `${r.genome_id} ${r.lineage} ${r.country}`.toLowerCase().includes(search));
  }
  // per-column value filters (the funnels in the header), same idea as the metadata table's quick filters
  const attrOf = (r) => { const g = meta.get(r.genome_id) || {};
    return { region: g.region || "", country: r.country || g.country || "", lineage: r.lineage || "", host: HOST_CLASS_LABEL[r.host_class] || "", bind: { ok: "정상", warn: "주의", bad: "미결합" }[bindState(r)] }; };
  state._msaAttrOf = attrOf; state._msaBase = rows;
  const colFilters = state.msaColFilters || {};
  rows = rows.filter((r) => Object.entries(colFilters).every(([k, set]) => set.has(attrOf(r)[k])));
  (state._msaRowsBySet = state._msaRowsBySet || new WeakMap()).set(setData, rows);   // all rows after search/filters -- the snapshot/copy buttons use these, not just the 800 drawn
  const shown = rows.slice(0, 800);

  const ampStartForRuler = setData.gene && setData.gene.amplicon_start !== "" && setData.gene.amplicon_start != null ? Number(setData.gene.amplicon_start) : NaN;
  const rulerPos = isNaN(ampStartForRuler) ? null : msaColumnToRefPos(setData, ampStartForRuler);
  const ruler = buildMsaRuler(msa.columns, rulerPos);
  const rulerTip = rulerPos
    ? "숫자 = 참조 유전체의 염기 좌표(10의 배수 위치마다 표시). 이 세트는 정렬이 좌표가 " + (rulerPos(0) > rulerPos(msa.columns - 1) ? "감소하는" : "증가하는") + " 방향으로 펼쳐져 있습니다. 삽입 열(참조에 없는 위치)에는 숫자가 없습니다."
    : "정렬의 열 번호 눈금 (10열마다)";

  const rowsHtml = shown.map((r) => {
    const mms = `<span class="msa-mm-tag ${mmClass(r.f_mismatch)}">F${mmText(r.f_mismatch)}</span>` +
      `<span class="msa-mm-tag ${mmClass(r.p_mismatch)}">P${mmText(r.p_mismatch)}</span>` +
      `<span class="msa-mm-tag ${mmClass(r.r_mismatch)}">R${mmText(r.r_mismatch)}</span>`;
    const g = meta.get(r.genome_id) || {};
    const bs = bindState(r);
    const bsLabel = { ok: "✅ 정상", warn: "⚠️ 주의", bad: "🚫 미결합" }[bs];
    const bsShort = { ok: "✅", warn: "⚠️", bad: "🚫" }[bs];
    const dash = (v) => (v === undefined || v === null || v === "" || v === "NA" ? "—" : escHtml(String(v)));
    return `<div class="msa-row msa-genome-row" data-genome-id="${r.genome_id}" data-label="${escHtml(`${r.genome_id} | ${r.lineage} | ${r.country || "?"}`)}" title="클릭하면 이 게놈의 결합 상세를 창으로 보여줍니다">
      ${metaCols
        ? `<div class="msa-label msa-label-cols"><span class="mono">${r.genome_id}</span><span>${dash(g.region)}</span><span>${dash(r.country || g.country)}</span><span>${dash(r.lineage)}</span><span>${dash(HOST_CLASS_LABEL[r.host_class])}</span>
        <span class="combo-status combo-status-${bs}">${bsLabel}</span><span>${mms}</span></div>`
        : `<div class="msa-label"><span class="mono">${bsShort} ${r.genome_id} | ${r.lineage} | ${r.country || "?"}</span> ${mms}</div>`}
      <div class="msa-seq">${markRanges(r.dots, ranges)}</div>
    </div>`;
  }).join("");

  const combined = buildCombinedSeqRow(msa);
  const strandParts = ["forward", "probe", "reverse"]
    .map((key) => {
      const rr = msa.ref_rows && msa.ref_rows[key];
      return rr ? `${key.toUpperCase()} ${strandLabel(rr.strand)}` : null;
    })
    .filter(Boolean)
    .join(" / ");
  const combinedRowHtml = `<div class="msa-row msa-ref-row msa-arrow-row msa-combined-row">
    <div class="msa-label"><b>PRIMER/PROBE 서열</b> <span class="msa-strand">(${strandParts})</span></div>
    <div class="msa-seq msa-arrow-combined">${markByClsArray(combined.text, combined.cls)}</div>
  </div>`;

  const excludedList = msa.excluded || [];
  const MSA_EXCLUDE_REASON_LABEL = {
    length_outlier: "길이 이상치 (amplicon 길이가 나머지 게놈들과 크게 다름)",
    disrupts_alignment: "정렬을 크게 왜곡시킴 (다른 게놈과 공유되지 않는 삽입/결손으로 전체 정렬을 망가뜨림)",
  };
  const excludedHtml = excludedList.length ? `
    <details class="msa-excluded-block">
      <summary>정렬에서 제외되고 별도로 나열된 게놈 ${excludedList.length}개 (문제가 있어 정렬 전체를 왜곡시키므로 분석에서 제외) &ndash; 펼쳐서 목록 보기</summary>
      <div class="msa-excluded-list">
        <div class="msa-excluded-row msa-excluded-head"><span>Accession</span><span>Lineage</span><span>제외 사유</span></div>
        ${excludedList.map((e) => `<div class="msa-excluded-row">
          <span class="mono">${e.genome_id}</span><span>${e.lineage || "NA"}</span>
          <span class="note">${MSA_EXCLUDE_REASON_LABEL[e.reason] || e.reason}</span>
        </div>`).join("")}
      </div>
    </details>` : "";

  return `
    <div class="table-toolbar">
      <input id="msaSearch" type="text" placeholder="accession, lineage, country 검색..." value="${state.msaSearch || ""}">
      ${Object.keys(state.msaColFilters || {}).length ? `<button type="button" id="msaClearColFilters" class="btn btn-secondary">열 필터 해제 (${Object.keys(state.msaColFilters).length})</button>` : ""}
      <label class="note" style="display:inline-flex;gap:4px;align-items:center;" title="대륙·나라·클레이드·결합 여부 열을 표시합니다 (좁은 화면에서는 끄면 정렬이 더 넓게 보입니다)"><input type="checkbox" id="msaMetaCols" ${metaCols ? "checked" : ""}> 정보 열</label>
      <select id="msaFilterMode" title="결합 판정은 게놈 표와 같은 기준입니다: F/R mismatch&le;1(3′ 끝 일치)·Probe 완전일치 = 결합, 그 가운데 Tm·&Delta;G·amplicon 크기가 의심스러우면 주의">
        <option value="all" ${mode === "all" ? "selected" : ""}>전체 (${msa.rows.length})</option>
        <option value="bound" ${mode === "bound" ? "selected" : ""}>결합 판정만: 정상 + 주의 (${counts.ok + counts.warn})</option>
        <option value="ok" ${mode === "ok" ? "selected" : ""}>정상만 (${counts.ok})</option>
        <option value="warn" ${mode === "warn" ? "selected" : ""}>주의만 (${counts.warn})</option>
        <option value="bad" ${mode === "bad" ? "selected" : ""}>미결합만 (${counts.bad})</option>
      </select>
      <span class="count">${shown.length} / ${rows.length}개 표시 (전체 ${msa.rows.length}개${excludedList.length ? `, ${excludedList.length}개 제외됨` : ""})</span>
    </div>
    <div id="msaColSelectReadout" class="msa-col-select-readout" hidden>
      <span>컬럼 선택: <b class="msa-col-select-range"></b></span>
      <button id="msaColSelectCopyBtn" class="btn btn-secondary">📋 선택 영역 복사</button>
      <button id="msaColSelectClearBtn" class="btn btn-secondary">선택 해제</button>
    </div>
    ${excludedHtml}
    <div class="msa-outer${metaCols ? " msa-wide-label" : ""}">
      <div class="msa-row msa-ruler">${metaCols
        ? `<div class="msa-label msa-label-cols msa-label-head"><span>Accession</span><span>대륙${funnel("region")}</span><span>나라${funnel("country")}</span><span>${escHtml((PROFILE.lineage && PROFILE.lineage.label) || "Lineage")}${funnel("lineage")}</span><span>숙주${funnel("host")}</span><span>결합${funnel("bind")}</span><span>F/P/R mismatch</span></div>`
        : `<div class="msa-label" title="${rulerTip}">ref position</div>`}<div class="msa-seq" title="${rulerTip}">${markRanges(ruler, ranges)}</div></div>
      <div class="msa-row msa-diversity-row" id="msaDiversityRow" title="리니지 평균 pixy &#960; (population 전체 다양성) -- 참조유전체 결합위치를 기준으로 한 근사 매핑입니다.">
        <div class="msa-label">&#960; (다양성)</div>
        <div class="msa-seq"><span class="note">불러오는 중...</span></div>
      </div>
      ${combinedRowHtml}
      <div class="msa-row msa-consensus"><div class="msa-label">CONSENSUS</div><div class="msa-seq">${markRanges(msa.consensus, ranges)}</div></div>
      ${rowsHtml}
    </div>
    <div class="legend" style="margin-top:8px;">
      <span><span class="dot" style="background:rgba(76,120,168,.4)"></span>Forward</span>
      <span><span class="dot" style="background:rgba(242,166,90,.5)"></span>Probe</span>
      <span><span class="dot" style="background:rgba(84,162,74,.4)"></span>Reverse</span>
      <span id="msaDiversityLegendNote"><span class="dot" style="background:rgba(84,162,74,.55)"></span>&#8594;<span class="dot" style="background:rgba(228,87,86,.55)"></span> &#960;(다양성) 행: 초록=보존(낮은 π) &#8594; 빨강=다양(높은 π), pixy 리니지 평균, 근사 매핑</span>
      <span>"." = consensus와 동일 (일치), 글자 = 미스매치, "-" = gap/indel</span>
      <span>PRIMER/PROBE 서열 행: <u style="text-decoration-style:solid;">밑줄</u> = 해당 위치가 population consensus와 불일치</span>
      <span>게놈 라벨: ACCESSION | LINEAGE | COUNTRY</span>
      <span>⌥/Alt + 드래그: VS Code처럼 특정 구간(컬럼)을 세로로 선택 (전체 행 비교용)</span>
    </div>
  `;
}

// the pass rule of the pipeline (Parameters page): F/R mismatch <= primer_mismatch_pass, probe mismatch <= probe_mismatch
const PASS = { mm: 1, pm: 0 };
function isValidatedGenome(g) {
  return g.f_mismatch !== null && g.f_mismatch <= PASS.mm &&
    g.p_mismatch !== null && g.p_mismatch <= PASS.pm &&
    g.r_mismatch !== null && g.r_mismatch <= PASS.mm;
}

// A lineage with very few sampled genomes in the whole filtered population gives coverage
// percentages and pixy pi values that look precise ("5/5, 100%") but are statistically
// close to meaningless -- flag it wherever a per-lineage stat is shown, rather than let a
// tiny sample look as authoritative as a well-sampled one.
const LOW_N_THRESHOLD = 20;
function lowNBadge(n) {
  if (n == null || n >= LOW_N_THRESHOLD) return "";
  return `<span class="pill pill-warn low-n-badge" title="population 내 이 리니지의 전체 게놈 수가 ${n}개로 적어(기준: ${LOW_N_THRESHOLD}개 미만) 이 값의 통계적 신뢰도가 낮습니다.">n=${n}</span>`;
}

function renderLineageBreakdownCompact(setData) {
  // Label+bar+value stacked vertically (rather than one wide row) -- this card sits in
  // the same card-grid row as the other 4 (narrower, ~200-260px) summary cards, where a
  // side-by-side label|bar|value layout has no room left for the bar itself.
  const genomes = setData.genomes || [];
  const totals = setData.lineage_totals || {};
  const counts = {};
  for (const g of genomes) {
    if (!isValidatedGenome(g)) continue;
    counts[g.lineage] = (counts[g.lineage] || 0) + 1;
  }
  // Every lineage present in the population (not just the ones this set happens to
  // validate at least one genome for) -- a lineage this set completely misses is exactly
  // the kind of gap worth seeing, not silently dropping from the list.
  const lineages = Array.from(new Set([...Object.keys(totals), ...Object.keys(counts)]));
  if (!lineages.length) return "<p class='note'>population 리니지 정보가 없습니다.</p>";
  const mode = state.lineageSortMode || "name";
  lineages.sort(mode === "name" ? (a, b) => a.localeCompare(b) : (a, b) => (counts[b] || 0) - (counts[a] || 0));
  return lineages.map((l) => {
    const denom = totals[l] || 0;
    const c = counts[l] || 0;
    const pct = denom ? (c / denom) * 100 : null;
    const barPct = denom ? Math.min(100, (c / denom) * 100) : 0;       // bar = the percentage printed beside it
    return `<div class="lineage-compact-row">
      <div class="lineage-compact-top"><span class="label">${l} ${lowNBadge(denom)}</span><span class="val">${c}${denom ? ` / ${denom} (${pct.toFixed(0)}%)` : ""}</span></div>
      <div class="cov-bar-track"><div class="cov-bar-fill validated" style="width:${barPct}%"></div></div>
    </div>`;
  }).join("");
}

const MAIN_VIEW_IDS = ["engineView", "overview", "detail", "addForm", "paramsView", "changelogView", "guideView", "diversityView", "geneConsView", "redesignView", "vennView", "metadataView", "primersView", "litView", "standardsView", "metaTableView", "stdMgrView", "requirementsView"];
function showMainView(id) {
  for (const v of MAIN_VIEW_IDS) $(`#${v}`).hidden = (v !== id);
  navRecord(id);
}

// ---- back / forward: every page change is a browser-history entry (Alt+arrows and the browser buttons work) and the bar on top of the page
// goes back / forward the same way. A page is {view id, active set/group id}; going back re-opens it through its own menu button.
const NAV = { idx: 0, max: 0, lock: false, last: null };
const NAV_VIEW_BTN = { overview: "overviewBtn", addForm: "addSetBtn", paramsView: "paramsBtn", changelogView: "changelogBtn", guideView: "guideBtn",
  diversityView: "diversityBtn", geneConsView: "geneConsBtn", redesignView: "redesignBtn", vennView: "vennBtn", metadataView: "metadataBtn",
  primersView: "primersBtn", litView: "litBtn", standardsView: "standardsBtn", metaTableView: "metaTableBtn", stdMgrView: "stdMgrBtn", requirementsView: "reportBtn" };
NAV_VIEW_BTN.requirementsView = "requirementsBtn";
NAV_VIEW_BTN.engineView = "engineBtn";
const NAV_TITLE = { engineView: "엔진 검증 (Go)", overview: "전체 대시보드", addForm: "새 프라이머 세트 분석", paramsView: "분석 파라미터", changelogView: "변경 이력", guideView: "파이프라인 개요 / 가이드",
  diversityView: "리니지별 유전체 다양성", geneConsView: "유전자 보존영역", redesignView: "리니지별 재설계", vennView: "커버리지 벤 다이어그램", metadataView: "메타데이터 / 지도",
  primersView: "프라이머 입력 / 관리", litView: "문헌 프라이머", standardsView: "표준물질 대표성", metaTableView: "메타데이터 표", stdMgrView: "표준물질 입력 / 관리", requirementsView: "요구사항 기록" };
function navTitle(d) {
  if (!d) return "";
  if (d.v === "detail") {
    if (d.at === "group") return "그룹 · " + ((state.groups.find((g) => g.group_id === d.aid) || {}).label || d.aid);
    return "세트 · " + ((state.sets.find((x) => x.set_id === d.aid) || {}).display_name || d.aid);
  }
  return NAV_TITLE[d.v] || "";
}
function navUpdate() {
  const b = document.getElementById("navBack"), f = document.getElementById("navFwd"), c = document.getElementById("navCrumb");
  if (!b) return;
  b.disabled = NAV.idx <= 0; f.disabled = NAV.idx >= NAV.max;
  c.textContent = navTitle(NAV.last);
}
function navRecord(id) {
  if (NAV.lock) return;
  const d = { v: id, aid: state.activeId || null, at: state.activeType || null };
  const same = NAV.last && NAV.last.v === d.v && NAV.last.aid === d.aid && NAV.last.at === d.at;
  if (same) { history.replaceState({ ...d, idx: NAV.idx }, ""); }
  else if (NAV.last) { NAV.idx += 1; NAV.max = NAV.idx; history.pushState({ ...d, idx: NAV.idx }, ""); }
  else { history.replaceState({ ...d, idx: NAV.idx }, ""); }
  NAV.last = d;
  navUpdate();
}
async function navRestore(st) {
  if (st.v === "detail") { return st.at === "group" ? selectGroup(st.aid) : selectSet(st.aid); }
  const btn = NAV_VIEW_BTN[st.v] && document.getElementById(NAV_VIEW_BTN[st.v]);
  if (btn) btn.click();
}
window.addEventListener("popstate", async (e) => {
  const st = e.state;
  if (!st || !st.v) return;
  NAV.lock = true;
  NAV.idx = st.idx || 0; NAV.last = { v: st.v, aid: st.aid, at: st.at };
  try { await navRestore(st); } catch (err) { /* page failed to open: stay on what is shown */ }
  await new Promise((r) => setTimeout(r, 150));
  NAV.lock = false; navUpdate();
});
document.addEventListener("DOMContentLoaded", () => {
  const prev = history.state;
  if (prev && prev.v) { NAV.idx = prev.idx || 0; NAV.max = NAV.idx; }
  document.getElementById("navBack").addEventListener("click", () => history.back());
  document.getElementById("navFwd").addEventListener("click", () => history.forward());
});

function narrateSetAssessment(setData) {
  // A heuristic, generic summary synthesized from the metrics already computed elsewhere
  // on this page (coverage, probe dropout, standard-material match, reference-genome
  // binding, amplicon size, MSA disruptive-genome count) -- meant as a quick "should I
  // be worried about this set" pointer, not a substitute for expert/wet-lab review.
  const s = setData.summary || {};
  const std = setData.standard || {};
  const gene = setData.gene || {};
  const binding = setData.binding || [];
  const msa = setData.msa || {};
  const points = [];
  let score = 0;

  const vpct = s.validated_pct ?? 0;
  if (vpct >= 95) { points.push({ ok: true, text: `Population 커버리지 ${fmtPct(vpct)} — 현재 알려진 유전적 다양성 대부분을 커버합니다.` }); score += 2; }
  else if (vpct >= 80) { points.push({ ok: true, text: `Population 커버리지 ${fmtPct(vpct)} — 양호한 편입니다.` }); score += 1; }
  else if (vpct >= 50) { points.push({ ok: null, text: `Population 커버리지 ${fmtPct(vpct)} — 절반 안팎에 그쳐, 상당수 lineage/변이체를 놓칠 수 있습니다.` }); score -= 1; }
  else { points.push({ ok: false, text: `Population 커버리지 ${fmtPct(vpct)} — 낮은 편이라 실제 시료의 상당수를 검출하지 못할 위험이 있습니다.` }); score -= 2; }

  const dropout = s.probe_dropout ?? 0;
  const primerMatched = s.primer_matched ?? 0;
  const dropoutFrac = primerMatched ? dropout / primerMatched : 0;
  if (dropoutFrac > 0.3) { points.push({ ok: false, text: `Primer는 매치되지만 Probe가 실패하는 게놈이 ${dropout}개(약 ${Math.round(dropoutFrac * 100)}%)로 많습니다 — Probe 서열(축퇴염기 보강 등) 재검토를 권장합니다.` }); score -= 1; }
  else if (dropoutFrac > 0.1) { points.push({ ok: null, text: `Probe dropout ${dropout}개(약 ${Math.round(dropoutFrac * 100)}%) — 일부 변이체에서 Probe 결합이 실패합니다.` }); }
  else if (primerMatched) { points.push({ ok: true, text: `Probe dropout ${dropout}개로 적어, Primer와 Probe 설계가 서로 잘 맞습니다.` }); score += 1; }

  if (std.available === false) { points.push({ ok: null, text: "표준물질이 등록되어 있지 않아 실험실 양성대조군과의 비교는 평가에서 제외했습니다(표준물질은 선택 사항)." }); }
  else if (std.ok) { points.push({ ok: true, text: `표준물질(KDCA) 매치 O (${(std.matched || []).join(", ")}) — 실제 wet-lab 양성대조군에서도 동작할 것으로 기대됩니다.` }); score += 1; }
  else { points.push({ ok: false, text: "표준물질(KDCA) 매치 X — 실험실 양성대조군 검증에 실패할 가능성이 있습니다." }); score -= 2; }

  const notFound = binding.filter((b) => !b.found);
  if (notFound.length) {
    points.push({ ok: false, text: `참조유전체에서 ${notFound.map((b) => b.role).join("/")} 결합위치를 찾지 못했습니다 — 서열 설계 자체를 재검토할 필요가 있습니다.` });
    score -= 2;
  } else if (gene.amplicon_size) {
    const size = Number(gene.amplicon_size);
    if (size > 300) { points.push({ ok: null, text: `Amplicon 크기 ${size}bp — TaqMan qPCR 권장 범위(약 60~150bp)보다 큰 편입니다.` }); }
    else { points.push({ ok: true, text: `Amplicon 크기 ${size}bp — qPCR에 적합한 범위입니다.` }); score += 1; }
  }

  const disruptiveCount = (msa.excluded || []).filter((e) => e.reason === "disrupts_alignment").length;
  if (disruptiveCount > 0) {
    points.push({ ok: null, text: `전체 게놈 정렬에서 ${disruptiveCount}개 게놈이 정렬을 크게 왜곡시켜 별도로 제외되었습니다 — 해당 게놈들은 개별적으로 결합 여부를 확인해보는 것이 좋습니다.` });
  }

  const dh = setData.dimer_hairpin;
  if (dh) {
    const threshold = dh.dg_flag_threshold ?? -9;
    const flaggedDimers = Object.entries(dh.dimer || {}).filter(([, v]) => v.worst_dg <= threshold);
    const flaggedHairpins = Object.entries(dh.hairpin || {}).filter(([, v]) => v.worst_dg <= threshold);
    if (flaggedDimers.length || flaggedHairpins.length) {
      const parts = [
        ...flaggedDimers.map(([pair, v]) => `${DIMER_PAIR_LABEL[pair]?.replace(/&harr;|&middot;/g, "-") || pair}(&#916;G ${v.worst_dg.toFixed(1)})`),
        ...flaggedHairpins.map(([role, v]) => `${DIMER_ROLE_LABEL[role] || role} hairpin(&#916;G ${v.worst_dg.toFixed(1)})`),
      ];
      points.push({ ok: false, text: `이차구조 위험: ${parts.join(", ")} — 자가/교차 결합이나 헤어핀 형성으로 증폭 효율이 떨어질 수 있습니다.` });
      score -= 1;
    } else {
      points.push({ ok: true, text: "프라이머/프로브 간 뚜렷한 self-dimer·cross-dimer·hairpin 위험은 발견되지 않았습니다." });
      score += 1;
    }
  }

  // Only a group's setData carries this -- an individual set has no "other member sets" to
  // cross-check against, so this key is simply absent there and the check is skipped, same
  // pattern as gene/binding above.
  const crossBinding = setData.cross_binding;
  if (crossBinding != null) {
    const crossGenomes = new Set(crossBinding.filter((g) => !g.is_same_set && g.genome_id).map((g) => g.genome_id));
    if (crossGenomes.size) {
      points.push({ ok: false, text: `그룹 구성 세트 간 의도치 않은 교차 결합이 population ${crossGenomes.size}개 게놈에서 발견됐습니다 — 비특이적 증폭(off-target) 가능성이 있으니 조합을 재검토하세요.` });
      score -= 2;
    } else {
      points.push({ ok: true, text: "그룹 구성 세트 간 의도치 않은 교차 결합은 population 전체에서 발견되지 않았습니다." });
      score += 1;
    }
  }

  let verdict, verdictClass;
  if (score >= 3) { verdict = "✅ 전반적으로 설계가 양호합니다."; verdictClass = "ok"; }
  else if (score >= 0) { verdict = "ℹ️ 대체로 사용 가능하나 일부 보완이 필요합니다."; verdictClass = "warn"; }
  else { verdict = "⚠️ 여러 문제가 발견되어 재설계 또는 보완을 권장합니다."; verdictClass = "bad"; }

  return { verdict, verdictClass, points };
}

// 요약(종합 평가) 오른쪽에 그리는 리니지별 누적 커버리지 곡선: 그룹은 구성 세트, 세트는 자기 자신 + 변이 세트를 하나씩 더할 때의 커버리지.
function summaryWithCurve(assessmentHtml, seg, keys) {
  if (!keys || keys.length < 2) return assessmentHtml;
  setTimeout(mountCurveCards, 0);
  return `<div class="summary-split">${assessmentHtml}<div class="card cov-curve-card" data-seg="${escHtml(String(seg).toLowerCase())}" data-keys="${escHtml(JSON.stringify(keys))}"><h3>리니지별 누적 커버리지 곡선</h3><p class="note">불러오는 중...</p></div></div>`;
}
async function mountCurveCards() {
  for (const el of document.querySelectorAll(".cov-curve-card:not([data-done])")) {
    el.dataset.done = "1";
    try {
      const keys = JSON.parse(el.dataset.keys);
      const d = await fetchJSON(`/api/venn/${el.dataset.seg}`);
      const byKey = new Map(d.items.map((it) => [it.key, it]));
      const sel = keys.map((k) => byKey.get(k)).filter(Boolean);
      if (sel.length < 2) { el.remove(); continue; }
      const data = vnCurveDataFrom(d.genomes, sel);
      const order = data.steps.slice(1).map((s) => `<li><b>${escHtml(s.label)}</b> <span class="note">+${s.gain}개 → 전체 ${(100 * s.cov.ALL / (data.total.ALL || 1)).toFixed(1)}%</span></li>`).join("");
      const body = `${vnCurveSvg(data)}<details style="margin-top:6px;"><summary class="note">추가 순서</summary><ol style="margin:6px 0 0 18px;">${order}</ol></details>`;
      el.innerHTML = `<h3 style="display:flex;justify-content:space-between;gap:6px;align-items:center;">리니지별 누적 커버리지 곡선 <button type="button" class="btn btn-secondary cov-curve-big" title="크게 보기">&#9974;</button></h3>
        <div class="note" style="margin-bottom:4px;">위로 빨리 올라갈수록 적은 수의 세트로 많이 잡습니다 (세트 단독 통과 기준)</div>${body}`;
      el.querySelector(".cov-curve-big").addEventListener("click", () => openDetailModal({ title: "리니지별 누적 커버리지 곡선", html: body }));
    } catch (e) { el.remove(); }
  }
}

function renderSetAssessment(setData) {
  const { verdict, verdictClass, points } = narrateSetAssessment(setData);
  return `
    <div class="card assessment-card ${verdictClass}" style="margin-bottom:20px;">
      <h3>종합 평가 (참고용 자동 분석 의견)</h3>
      <div class="assessment-verdict ${verdictClass}">${verdict}</div>
      <ul class="assessment-points">
        ${points.map((p) => `<li class="${p.ok === true ? "ok" : p.ok === false ? "bad" : "warn"}">${p.ok === true ? "✅" : p.ok === false ? "⚠️" : "ℹ️"} ${p.text}</li>`).join("")}
      </ul>
      <div class="note" style="margin-top:8px;">이 평가는 계산된 지표를 바탕으로 한 자동 요약이며, 최종 판단은 실험 검증과 전문가 리뷰를 거쳐야 합니다.</div>
    </div>`;
}

function renderDetail(setData) {
  state.detailData = setData;      // read by the standard-material pop-up (#stdMatchCard)
  showMainView("detail");
  const detail = $("#detail");

  const s = setData.summary || {};
  const std = setData.standard || {};
  const gene = setData.gene || {};

  const setDisplayName = setData.display_name || setData.set_id;
  detail.innerHTML = `
    <div class="detail-header">
      <h2>${escHtml(setDisplayName)}</h2>
      <button id="renameSetBtn" class="btn btn-secondary" title="이름 변경">&#9998; 이름 변경</button>
      ${setData.builtin === false ? "" : `<button id="analyzeSetBtn" class="btn btn-secondary" title="이 세트만 분석합니다. 이미 최신인 단계는 건너뛰므로 새 세트나 바뀐 세트에 빠릅니다.">&#9654; 이 세트 분석</button>
      <button id="reanalyzeSetBtn" class="btn btn-secondary" title="입력이 같아도 이 세트의 결과를 처음부터 다시 계산합니다 (다른 세트는 건드리지 않음).">&#8635; 이 세트 재분석</button>`}
      <span class="pill pill-${setData.segment}">${setData.segment}</span>
      <span class="pill ${setData.builtin === false ? "pill-user" : ""}">${setData.source}</span>
    </div>
    ${finalResultBoxHtml(setData.set_id)}
    ${relationsHtml(setData.set_id)}
    <div class="detail-sub">
      ${setDisplayName !== setData.set_id ? `ID: <span class="mono">${escHtml(setData.set_id)}</span> &middot; ` : ""}Target gene: <b>${gene.target_gene || "NA"}</b>
      ${gene.amplicon_start ? ` &middot; amplicon ${gene.chrom || ""}:${gene.amplicon_start}..${gene.amplicon_end} (${gene.amplicon_size}bp)` : ""}
    </div>

    ${summaryWithCurve(renderSetAssessment(setData), setData.segment, (() => { const vs = state.sets.filter((x) => x.parent_id === setData.set_id).map((x) => x.set_id); return vs.length ? [setData.set_id, ...vs] : []; })())}

    <div class="card-grid">
      <div class="card">
        <h3>프라이머 / 프로브 서열</h3>
        <div class="seq-row"><span class="tag">Forward</span><span class="seq">${setData.F_seq}</span></div>
        <div class="seq-row"><span class="tag">Probe</span><span class="seq">${setData.P_seq}</span></div>
        <div class="seq-row"><span class="tag">Reverse</span><span class="seq">${setData.R_seq}</span></div>
        <div class="seq-row" style="margin-top:8px;align-items:center;">
          <span class="tag">형광물질</span>
          <select id="fluorophoreSelect" style="padding:4px 6px;border-radius:6px;border:1px solid var(--border);background:var(--bg);color:var(--fg);font-size:12px;">
            ${renderFluorophoreOptions(setData.fluorophore)}
          </select>
          <input id="fluorophoreCustom" type="text" placeholder="형광물질 이름 직접 입력" value="${!FLUOROPHORE_PRESETS.includes(setData.fluorophore) ? escHtml(setData.fluorophore || "") : ""}" ${!FLUOROPHORE_PRESETS.includes(setData.fluorophore) && setData.fluorophore ? "" : "hidden"} style="padding:4px 6px;border-radius:6px;border:1px solid var(--border);background:var(--bg);color:var(--fg);font-size:12px;width:120px;">
          <button id="saveFluorophoreBtn" class="btn btn-secondary" style="padding:4px 10px;font-size:11px;">저장</button>
          <span id="fluorophoreStatus" class="note"></span>
        </div>
        <div class="note" style="margin-top:4px;">실험(Wet-lab) 검증 시 어느 채널(형광 색상)의 증폭곡선과 비교해야 하는지 알려줍니다.</div>
        ${renderCitation(setData.citation)}
      </div>
      <div class="card">
        <h3>Population 커버리지</h3>
        <div class="metric-list">
          <div class="row"><span>Primer mismatch=0</span><span class="v">${fmtPct(s.primer_coverage_mm0_pct)}</span></div>
          <div class="row"><span>Primer mismatch&le;1</span><span class="v">${fmtPct(s.primer_coverage_mm1_pct)}</span></div>
          <div class="row"><span>Primer mismatch&le;2</span><span class="v">${fmtPct(s.primer_coverage_mm2_pct)}</span></div>
          <div class="row"><span>Probe 완전일치</span><span class="v">${fmtPct(s.probe_coverage_mm0_pct)}</span></div>
          <div class="row"><span>Probe dropout</span><span class="v">${s.probe_dropout ?? "—"}</span></div>
          <div class="row" style="border-top:1px solid var(--border);padding-top:6px;margin-top:2px;">
            <span><b>최종 Validated</b></span><span class="v" style="color:var(--accent);">${s.validated ?? "—"} / ${s.total_genomes ?? "—"} (${fmtPct(s.validated_pct)})</span></div>
        </div>
      </div>
      <div class="card clickable-card" id="stdMatchCard" role="button" tabindex="0" title="클릭하면 표준물질 매치 상세를 창으로 보여줍니다">
        <h3>표준물질(KDCA) 매치</h3>
        <div class="big">${std.available === false ? "<span style='color:var(--muted)'>&mdash;</span>" : std.ok ? "<span style='color:var(--ok)'>O</span>" : "<span style='color:var(--bad)'>X</span>"}
          <small>${std.ok ? (std.matched || []).join(", ") : ""}</small></div>
        ${std.available === false ? `<div class="note" style="margin-top:6px;">표준물질 미등록(선택 사항)</div>` : std.ok ? `<div class="note" style="margin-top:6px;">${std.perfect_primer ? "완전일치" : `primer mismatch&le;1 허용`}</div>` : ""}
      </div>
      ${renderStdRepCard(setData)}
      <div class="card clickable-card" id="bindingCard" role="button" tabindex="0" title="클릭하면 참조유전체 결합 상세를 창으로 보여줍니다">
        <h3>참조유전체 결합위치</h3>
        ${renderBinding(setData.binding)}
      </div>
      <div class="card">
        <h3>이차구조 (Self/Cross-dimer &middot; Hairpin)</h3>
        ${renderDimerHairpin(setData.dimer_hairpin)}
      </div>
      <div class="card">
        <h3 style="display:flex;align-items:center;justify-content:space-between;gap:6px;">
          검증 통과 게놈의 리니지 분포
          <button id="lineageSortBtn" class="sort-toggle-btn" title="정렬 방식 전환">⇅ <span id="lineageSortLabel">${state.lineageSortMode === "name" ? "이름순" : "건수순"}</span></button>
        </h3>
        <div id="lineageCardBody">${renderLineageBreakdownCompact(setData)}</div>
      </div>
    </div>

    ${renderSetDiversitySection(setData)}

    <div class="card" style="margin-bottom:20px;">
      <h3>표준물질(KDCA) 상세 결합 정보</h3>
      <div class="note" style="margin-bottom:10px;">
        표준물질 각각에 대해 F/Probe/R이 정확히 어디에 어떻게 결합하는지, 매치가 안 되는 경우
        mismatch가 몇 개인지까지 게놈별 미스매치 표와 동일한 방식으로 표시합니다.
      </div>
      ${renderStandardDetail(setData.standard_genomes, setData)}
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3 style="display:flex;align-items:center;justify-content:space-between;gap:6px;">
        🧪 실험(Wet-lab) 검증
        <button id="openVerifyModalBtn" class="btn btn-primary">+ 검증 결과 추가</button>
      </h3>
      <div class="note" style="margin-bottom:10px;">
        표준물질(양성대조물질)에 실제 프라이머/프로브를 넣어 테스트한 RT-PCR 결과(증폭곡선, 겔 사진 등)를
        업로드하고, 위 인실리코 예측(mismatch, Tm, ΔG, O/X 매치)과 비교해 최종 검증 상태를 기록합니다.
      </div>
      <div id="verifyHistoryOnPage"></div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>메모 / 노트</h3>
      <textarea id="setNoteInput" rows="3" placeholder="이 프라이머 세트에 대한 메모를 남겨두세요 (예: 설계 근거, 검토 이력, 주의사항)... 사진은 여기에 붙여넣기(Ctrl/Cmd+V)하거나 아래 버튼으로 추가할 수 있습니다." style="width:100%;resize:vertical;padding:8px 10px;border-radius:8px;border:1px solid var(--border);background:var(--bg);color:var(--fg);font-size:12.5px;">${setData.note || ""}</textarea>
      <div style="display:flex;align-items:center;gap:10px;margin-top:6px;">
        <button id="saveNoteBtn" class="btn btn-primary">저장</button>
        <button id="addNoteImageBtn" class="btn btn-secondary">📷 사진 추가 (예: RT-PCR 결과)</button>
        <input id="noteImageFileInput" type="file" accept="image/png,image/jpeg,image/gif,image/webp" multiple hidden>
        <span id="noteSaveStatus" class="note"></span>
      </div>
      <div id="noteImageGallery" class="note-image-gallery"></div>
    </div>

    <div class="card msa-card" style="margin-bottom:20px;">
      <div class="msa-card-header">
        <h3 style="margin:0;">전체 게놈 정렬 (F / Probe / R 결합을 한번에 보기)</h3>
        <div style="display:flex;gap:6px;flex-wrap:wrap;"><button id="msaSnapBtn" class="btn btn-secondary btn-fullscreen" title="현재 필터가 적용된 정렬 전체(잘리는 부분 없이)를 PNG 이미지로 저장합니다. 컬럼을 선택해 두었으면 그 구간만 저장합니다.">&#128247; 이미지 저장</button>
        <button id="msaCopyImgBtn" class="btn btn-secondary btn-fullscreen" title="같은 이미지를 클립보드로 복사합니다 (HTTPS 또는 localhost에서만 가능, 아니면 저장됩니다)">&#128203; 이미지 복사</button>
        <button id="msaCopyTxtBtn" class="btn btn-secondary btn-fullscreen" title="정렬을 텍스트(accession | lineage | country 탭 서열)로 클립보드에 복사합니다. 컬럼을 선택해 두었으면 그 구간만 복사합니다.">&#128203; 텍스트 복사</button>
        <button id="msaFullscreenBtn" class="btn btn-secondary btn-fullscreen" title="전체화면">&#9974; 전체화면</button></div>
      </div>
      ${renderMSA(setData)}
    </div>

    <div class="card probe-design-card" id="probeDesignCard" style="margin-bottom:20px;">
      <div class="msa-card-header">
        <h3 style="margin:0;">Probe 자동 설계 · 테스트 <span class="note">(이 세트의 F/R 사이에서 후보 생성)</span></h3>
        <div style="display:flex;gap:6px;align-items:center;flex-wrap:wrap;">
          <label class="note" title="표시할 후보 수 (최적 탐색에서는 전부 검증합니다)">후보 <input id="probeDesignTop" type="number" min="1" max="30" value="20" style="width:56px;"></label>
          <label class="note" title="후보를 어디서 뽑을지: 자체 규칙 / Primer3의 internal oligo 설계 / 둘 다 (둘 다이면 절반씩)">후보 출처 <select id="probeDesignMethod"><option value="both">둘 다</option><option value="primer3">Primer3</option><option value="own">자체 규칙</option></select></label>
          <label class="note" title="Probe Tm을 두 프라이머 중 더 높은 쪽보다 몇 도 높게 할지 (최소 ~ 최대). 기본 +3 ~ +10 °C">Tm + <input id="probeDesignTmMin" type="number" step="0.5" value="3" style="width:52px;">~<input id="probeDesignTmMax" type="number" step="0.5" value="10" style="width:52px;">&deg;C</label>
          <label class="note" title="후보 1개에 넣을 수 있는 IUPAC(변이 허용, 예: R=A/G) 염기 수. 0이면 consensus 그대로">IUPAC 최대 <input id="probeDesignIupac" type="number" min="0" max="4" value="2" style="width:48px;"></label>
          <label class="note" title="MFEprimer population 검색의 최대 mismatch (기본 0 = 파이프라인의 Probe 기준)">mismatch <input id="probeDesignMm" type="number" min="0" max="3" value="0" style="width:48px;"></label>
          <button id="probeDesignRunBtn" class="btn btn-secondary">&#9654; Probe 설계 실행</button>
          <button id="probeBestBtn" class="btn btn-primary" title="위 조건으로 후보를 만들고, 후보 전부를 실제 검증 과정에 태워 결과가 가장 좋은 Probe를 찾습니다 (후보당 약 7초, 후보 20개면 2~3분). 후보 수는 위의 '후보' 칸을 따릅니다 (이 버튼은 기본 20).">&#127942; 최적 Probe 자동 탐색</button>
        </div>
      </div>
      <div id="probeDesignBody"><p class="note">불러오는 중...</p></div>
      <div id="probeDesignBest"></div>
      <div id="probeDesignValidations"></div>
    </div>

    <div class="card tree-card" style="margin-bottom:20px;">
      <div class="msa-card-header">
        <h3 style="margin:0;">Nextstrain 계통수 대비 커버리지</h3>
        <button id="treeFullscreenBtn" class="btn btn-secondary btn-fullscreen" title="전체화면">&#9974; 전체화면</button>
      </div>
      <div id="treeContainer">
        <button id="loadTreeBtn" class="btn btn-secondary">계통수 불러오기</button>
        <div class="note" style="margin-top:6px;">auspice.json(nextstrain.org/lassa/${setData.segment.toLowerCase()})의 원본 트리에서
        이 세트가 어떤 게놈까지 검증 가능한지 색으로 표시합니다.</div>
      </div>
    </div>

    <div class="card">
      <h3>게놈별 F / Probe / R 미스매치</h3>
      <div class="table-toolbar">
        <input id="genomeSearch" type="text" placeholder="accession, lineage, country 검색...">
        <label><input id="problemsOnlyChk" type="checkbox"> 문제/주의만 보기</label>
        <button id="genomeFilterResetBtn" class="btn btn-secondary" title="검색어, 필터, 정렬을 모두 초기 상태로">필터 리셋</button>
        <span id="genomeCount" class="count"></span>
      </div>
      <div class="genome-table-wrap">
        <table class="genome-table">
          <thead><tr>
            <th>상태</th>
            <th data-key="genome_id">Accession</th>
            <th data-key="lineage">Lineage</th>
            <th data-key="country">Country</th>
            <th data-key="f_mismatch">F mm</th>
            <th data-key="p_mismatch">P mm</th>
            <th data-key="r_mismatch">R mm</th>
          </tr></thead>
          <tbody id="genomeTableBody"></tbody>
        </table>
      </div>
      <div class="legend">
        <span><span class="combo-status combo-status-ok">✅ 정상</span></span>
        <span><span class="combo-status combo-status-warn">⚠️ 주의</span></span>
        <span><span class="combo-status combo-status-bad">🚫 문제</span></span>
        <span><span class="dot" style="background:rgba(84,162,74,.5)"></span>0 mismatch</span>
        <span><span class="dot" style="background:rgba(242,166,90,.6)"></span>1 mismatch</span>
        <span><span class="dot" style="background:rgba(228,87,86,.6)"></span>2+ mismatch</span>
        <span><span class="dot" style="background:var(--muted)"></span>no data</span>
        <span>행 클릭 시 amplicon 서열 (F=파랑, Probe=주황, R=초록) 표시</span>
      </div>
    </div>
  `;

  const genomes = setData.genomes || [];
  renderGenomeTable(genomes, setData);
  $("#genomeSearch").addEventListener("input", () => renderGenomeTable(genomes, setData));
  $("#problemsOnlyChk").addEventListener("change", (e) => { state.problemsOnly = e.target.checked; renderGenomeTable(genomes, setData); });
  $("#genomeFilterResetBtn").addEventListener("click", () => {
    $("#genomeSearch").value = "";
    $("#problemsOnlyChk").checked = false;
    state.problemsOnly = false;
    state.genomeSort = { key: "genome_id", dir: 1 };
    renderGenomeTable(genomes, setData);
  });
  detail.querySelectorAll("th[data-key]").forEach((th) => {
    th.addEventListener("click", () => {
      const key = th.dataset.key;
      state.genomeSort = { key, dir: state.genomeSort.key === key ? -state.genomeSort.dir : 1 };
      renderGenomeTable(genomes, setData);
    });
  });

  initProbeDesign(detail, setData);
  const msaCard = detail.querySelector(".msa-card");
  const MSA_HEADER = `<div class="msa-card-header">
      <h3 style="margin:0;">전체 게놈 정렬 (F / Probe / R 결합을 한번에 보기)</h3>
      <div style="display:flex;gap:6px;">
        <button id="msaExportBtn" class="btn btn-secondary btn-fullscreen" title="FASTA로 내보내서 SeaLion(전문 정렬 뷰어)에서 열기">&#128190; SeaLion용 FASTA</button>
        <button id="msaSnapBtn" class="btn btn-secondary btn-fullscreen" title="현재 필터가 적용된 정렬 전체(잘리는 부분 없이)를 PNG 이미지로 저장합니다. 컬럼을 선택해 두었으면 그 구간만 저장합니다.">&#128247; 이미지 저장</button>
        <button id="msaCopyImgBtn" class="btn btn-secondary btn-fullscreen" title="같은 이미지를 클립보드로 복사합니다 (HTTPS 또는 localhost에서만 가능, 아니면 저장됩니다)">&#128203; 이미지 복사</button>
        <button id="msaCopyTxtBtn" class="btn btn-secondary btn-fullscreen" title="정렬을 텍스트(accession | lineage | country 탭 서열)로 클립보드에 복사합니다. 컬럼을 선택해 두었으면 그 구간만 복사합니다.">&#128203; 텍스트 복사</button>
        <button id="msaFullscreenBtn" class="btn btn-secondary btn-fullscreen" title="전체화면">&#9974; 전체화면</button>
      </div>
    </div>`;
  const refreshMsa = () => {
    const prevFocus = document.activeElement && document.activeElement.id === "msaSearch";
    const caret = prevFocus ? document.activeElement.selectionStart : null;
    msaCard.innerHTML = MSA_HEADER + renderMSA(setData);
    wireMsa();
    // renderMSA() writes an empty pi row; refill it from the already-loaded diversity data (the first load fills it itself)
    if (diversityState.data) updateMsaDiversityBaseline(setData, diversityState.data);
    if (prevFocus) {
      const el2 = $("#msaSearch");
      el2.focus();
      el2.setSelectionRange(caret, caret);
    }
  };
  function wireMsa() {
    const searchInput = $("#msaSearch");
    if (searchInput) searchInput.addEventListener("input", (e) => { state.msaSearch = e.target.value; refreshMsa(); });
    msaCard.querySelectorAll(".mt-funnel[data-msa-col]").forEach((b) =>
      b.addEventListener("click", (e) => { e.stopPropagation(); msaOpenColMenu(b, b.dataset.msaCol, refreshMsa); }));
    const clearCols = $("#msaClearColFilters");
    if (clearCols) clearCols.addEventListener("click", () => { state.msaColFilters = {}; refreshMsa(); });
    const metaSel = $("#msaMetaCols");
    if (metaSel) metaSel.addEventListener("change", (e) => { state.msaMetaCols = e.target.checked; refreshMsa(); });
    const modeSel = $("#msaFilterMode");
    if (modeSel) modeSel.addEventListener("change", (e) => { state.msaFilterMode = e.target.value; refreshMsa(); });
    const fsBtn = $("#msaFullscreenBtn");
    if (fsBtn) fsBtn.addEventListener("click", () => toggleFullscreen(msaCard));
    const snapBtn = $("#msaSnapBtn");
    if (snapBtn) snapBtn.addEventListener("click", () => msaSnapshot(msaCard, setData, "save"));
    const copyImgBtn = $("#msaCopyImgBtn");
    if (copyImgBtn) copyImgBtn.addEventListener("click", () => msaSnapshot(msaCard, setData, "copy"));
    const copyTxtBtn = $("#msaCopyTxtBtn");
    if (copyTxtBtn) copyTxtBtn.addEventListener("click", () => msaCopyText(msaCard, setData));
    const exportBtn = $("#msaExportBtn");
    if (exportBtn) exportBtn.addEventListener("click", () => exportMsaFasta(setData));
  }
  wireMsa();
  // Clicking a genome in the MSA jumps to (and expands) that same genome's row in the
  // "게놈별 F/Probe/R 미스매치" table below -- attached once on msaCard itself (which
  // refreshMsa() never replaces, only its innerHTML), delegated to the re-rendered rows so
  // it isn't lost/duplicated across search or filter changes.
  msaCard.addEventListener("click", (e) => {
    const row = e.target.closest(".msa-genome-row");
    if (!row) return;
    const g = (setData.genomes || []).find((x) => x.genome_id === row.dataset.genomeId);
    if (g) openGenomeModal(g, setData);
  });
  wireColumnSelection(msaCard);
  document.addEventListener("fullscreenchange", () => {
    const fsBtn = $("#msaFullscreenBtn");
    if (fsBtn) fsBtn.innerHTML = document.fullscreenElement ? "&#9974; 전체화면 종료" : "&#9974; 전체화면";
    const tfsBtn = $("#treeFullscreenBtn");
    if (tfsBtn) tfsBtn.innerHTML = document.fullscreenElement ? "&#9974; 전체화면 종료" : "&#9974; 전체화면";
  });

  const treeCard = detail.querySelector(".tree-card");
  $("#loadTreeBtn").addEventListener("click", () => renderTreeSection($("#treeContainer"), setData, setData.segment));
  $("#treeFullscreenBtn").addEventListener("click", () => toggleFullscreen(treeCard));

  $("#lineageSortBtn").addEventListener("click", () => {
    state.lineageSortMode = state.lineageSortMode === "name" ? "count" : "name";
    $("#lineageSortLabel").textContent = state.lineageSortMode === "name" ? "이름순" : "건수순";
    $("#lineageCardBody").innerHTML = renderLineageBreakdownCompact(setData);
  });

  $("#saveNoteBtn").addEventListener("click", async () => {
    const statusEl = $("#noteSaveStatus");
    statusEl.textContent = "저장 중...";
    try {
      const note = $("#setNoteInput").value;
      const res = await fetch(`/api/sets/${encodeURIComponent(setData.set_id)}/note`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ note }),
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      const s = state.sets.find((x) => x.set_id === setData.set_id);
      if (s) s.note = note;
      statusEl.textContent = "저장됨";
      renderSidebar();
    } catch (e) {
      statusEl.textContent = `저장 실패: ${e.message}`;
    }
  });

  wireNoteImages(setData);
  wireVerification(setData);
  $("#renameSetBtn").addEventListener("click", () => renameSet(setData.set_id, setDisplayName));
  wireFluorophore(setData);
  loadSetDiversityOverlay(setData);
}

function wireFluorophore(setData) {
  const sel = $("#fluorophoreSelect");
  const custom = $("#fluorophoreCustom");
  sel.addEventListener("change", () => { custom.hidden = sel.value !== "__custom__"; if (!custom.hidden) custom.focus(); });
  $("#saveFluorophoreBtn").addEventListener("click", async () => {
    const statusEl = $("#fluorophoreStatus");
    const value = sel.value === "__custom__" ? custom.value.trim() : sel.value;
    statusEl.textContent = "저장 중...";
    try {
      const res = await fetch(`/api/sets/${encodeURIComponent(setData.set_id)}/fluorophore`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ fluorophore: value }),
      });
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
      setData.fluorophore = data.fluorophore;
      const s = state.sets.find((x) => x.set_id === setData.set_id);
      if (s) s.fluorophore = data.fluorophore;
      statusEl.textContent = "저장됨";
    } catch (e) {
      statusEl.textContent = `저장 실패: ${e.message}`;
    }
  });
}

// ---------------- note images (e.g. RT-PCR gel photos attached to a set's memo) ----------------

function renderNoteImageGallery(images, entityId, kind = "sets") {
  if (!images || !images.length) return "";
  return images.map((img) => `
    <div class="note-image-thumb" data-image-id="${img.id}">
      <img src="/api/${kind}/${encodeURIComponent(entityId)}/note/images/${img.id}" alt="${img.caption || "첨부 이미지"}" loading="lazy">
      <button class="note-image-delete" data-image-id="${img.id}" title="이미지 삭제">&times;</button>
      ${img.caption ? `<div class="note-image-caption">${img.caption}</div>` : ""}
    </div>`).join("");
}

function openImageLightbox(src) {
  const overlay = el("div", "modal-overlay lightbox-overlay");
  overlay.innerHTML = `<img src="${src}" class="lightbox-img">`;
  overlay.addEventListener("click", () => overlay.remove());
  modalHost().appendChild(overlay);
}

function wireNoteImages(setData, kind = "sets") {
  const gallery = $("#noteImageGallery");
  const setId = setData.set_id || setData.group_id;
  let images = (setData.note_images || []).slice();

  const refreshGallery = () => {
    gallery.innerHTML = renderNoteImageGallery(images, setId, kind);
    gallery.querySelectorAll(".note-image-thumb img").forEach((imgEl) => {
      imgEl.addEventListener("click", () => openImageLightbox(imgEl.src));
    });
    gallery.querySelectorAll(".note-image-delete").forEach((btn) => {
      btn.addEventListener("click", async (e) => {
        e.stopPropagation();
        const imageId = btn.dataset.imageId;
        if (!await uiConfirm("이 이미지를 삭제할까요?")) return;
        try {
          const res = await fetch(`/api/${kind}/${encodeURIComponent(setId)}/note/images/${imageId}`, { method: "DELETE" });
          if (!res.ok) throw new Error(`HTTP ${res.status}`);
          images = images.filter((img) => img.id !== imageId);
          refreshGallery();
        } catch (err) {
          uiAlert(`삭제 실패: ${err.message}`);
        }
      });
    });
  };
  refreshGallery();

  const uploadFile = async (file) => {
    if (!file.type.startsWith("image/")) return;
    const statusEl = $("#noteSaveStatus");
    statusEl.textContent = "이미지 업로드 중...";
    try {
      const dataUrl = await new Promise((resolve, reject) => {
        const reader = new FileReader();
        reader.onload = () => resolve(reader.result);
        reader.onerror = reject;
        reader.readAsDataURL(file);
      });
      const base64 = dataUrl.split(",")[1];
      const res = await fetch(`/api/${kind}/${encodeURIComponent(setId)}/note/images`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ content_type: file.type, data_base64: base64, caption: file.name }),
      });
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
      images.push(data.image);
      refreshGallery();
      statusEl.textContent = "이미지 추가됨";
    } catch (err) {
      $("#noteSaveStatus").textContent = `이미지 업로드 실패: ${err.message}`;
    }
  };

  $("#addNoteImageBtn").addEventListener("click", () => $("#noteImageFileInput").click());
  $("#noteImageFileInput").addEventListener("change", (e) => {
    [...e.target.files].forEach(uploadFile);
    e.target.value = "";
  });
  $("#setNoteInput").addEventListener("paste", (e) => {
    const items = [...(e.clipboardData?.items || [])];
    const imageItem = items.find((it) => it.type.startsWith("image/"));
    if (!imageItem) return;
    e.preventDefault();
    uploadFile(imageItem.getAsFile());
  });
}

// ---------------- wet-lab verification (real RT-PCR result vs in-silico prediction) ----------------
// A standard-material (양성대조물질) RT-PCR result the user ran at the bench, attached with photos
// (amplification curve, gel image, etc.) and compared against the already-computed in-silico
// prediction for that same material (setData.standard_genomes' mismatch/Tm/ΔG, via verdictFor()) --
// flags a discrepancy when the wet-lab call and the in-silico call disagree.

function verifyInSilicoFor(setData, materialId) {
  const g = (setData.standard_genomes || []).find((x) => x.genome_id === materialId);
  if (!g) return { g: null, ok: null };
  return { g, ok: verdictFor(g).ok };
}

function renderVerifyInSilicoPreview(setData, materialId) {
  if (!materialId) return "<p class='note'>표준물질을 선택하면 인실리코 예측(mismatch/Tm/ΔG)이 여기 표시됩니다.</p>";
  const { g, ok } = verifyInSilicoFor(setData, materialId);
  if (!g) return "<p class='note'>이 표준물질에 대한 인실리코 결과가 없습니다.</p>";
  return `
    <div class="verify-insilico-box">
      <div class="verify-insilico-title">인실리코 예측
        <span class="pill ${ok ? "pill-ok" : "pill-bad"}">${ok ? "O 매치" : "X 불일치"}</span>
      </div>
      ${highlightAmplicon(g, setData)}
    </div>`;
}

function verifyStatusPill(entry, setData) {
  const { ok: inSilicoOk } = verifyInSilicoFor(setData, entry.standard_material);
  const wetOk = entry.detected;
  const wetLabel = wetOk === true ? "실험: 검출(양성)" : wetOk === false ? "실험: 미검출(음성)" : "실험: 판정보류";
  const wetCls = wetOk === true ? "pill-ok" : wetOk === false ? "pill-bad" : "";
  const inSilicoLabel = inSilicoOk === true ? "인실리코: O 매치" : inSilicoOk === false ? "인실리코: X 불일치" : "인실리코: 데이터없음";
  const inSilicoCls = inSilicoOk === true ? "pill-ok" : inSilicoOk === false ? "pill-bad" : "";
  let compareHtml = "";
  if (inSilicoOk !== null && wetOk !== null) {
    compareHtml = inSilicoOk === wetOk
      ? `<span class="pill pill-ok">✅ 인실리코-실험 일치</span>`
      : `<span class="pill pill-warn">⚠ 인실리코-실험 불일치 — 확인 필요</span>`;
  }
  return `<span class="pill ${wetCls}">${wetLabel}</span><span class="pill ${inSilicoCls}">${inSilicoLabel}</span>${compareHtml}`;
}

function renderVerificationHistory(entries, setId) {
  if (!entries || !entries.length) return "<p class='note'>아직 등록된 실험 검증 기록이 없습니다.</p>";
  return entries.slice().reverse().map((entry) => `
    <div class="verify-entry" data-vid="${entry.id}">
      <div class="verify-entry-head">
        <span class="mono">${escHtml(entry.standard_material)}</span>
        ${entry.fluorophore ? `<span class="pill" title="형광물질(채널)">${escHtml(entry.fluorophore)}</span>` : ""}
        <span id="vhPills-${entry.id}"></span>
        <span class="note">${escHtml(entry.created_at)}${entry.ct_value ? ` &middot; Ct ${escHtml(entry.ct_value)}` : ""}</span>
        <button class="btn btn-secondary verify-delete-btn" data-vid="${entry.id}" title="검증 기록 삭제">삭제</button>
      </div>
      ${entry.note ? `<div class="note" style="margin:4px 0;">${escHtml(entry.note)}</div>` : ""}
      <div class="note-image-gallery" data-vid="${entry.id}" id="vhGallery-${entry.id}"></div>
      <div style="margin-top:6px;">
        <button class="btn btn-secondary verify-add-image-btn" data-vid="${entry.id}">📷 결과 이미지 추가</button>
        <input type="file" class="verify-image-file-input" data-vid="${entry.id}" accept="image/png,image/jpeg,image/gif,image/webp" multiple hidden>
      </div>
    </div>`).join("");
}

// #vfSaveBtn and #vfMaterial live inside the page-level #verifyModal singleton (not recreated
// per render), so their listeners must be registered once at module load and read whichever
// set is currently open via this shared context, instead of being re-registered (and piling up)
// on every wireVerification() call. See the module-level listeners below wireVerification.
let currentVerifyTarget = null;

function wireVerification(setData) {
  const setId = setData.set_id;
  let entries = (setData.verifications || []).slice();

  const materialOptions = () => (setData.standard_genomes || []).map((g) => {
    const ok = verdictFor(g).ok;
    return `<option value="${g.genome_id}">${g.genome_id} (인실리코: ${ok ? "O" : "X"})</option>`;
  }).join("");

  function refreshOnPage() {
    const target = $("#verifyHistoryOnPage");
    if (!target) return;
    target.innerHTML = renderVerificationHistory(entries, setId);
    wireHistoryList(target);
  }

  function wireHistoryList(container) {
    entries.forEach((entry) => {
      const pillsEl = container.querySelector(`#vhPills-${entry.id}`);
      if (pillsEl) pillsEl.innerHTML = verifyStatusPill(entry, setData);
      const galleryEl = container.querySelector(`#vhGallery-${entry.id}`);
      if (galleryEl) {
        galleryEl.innerHTML = renderNoteImageGalleryFor(entry.images, setId, entry.id);
        galleryEl.querySelectorAll(".note-image-thumb img").forEach((imgEl) => {
          imgEl.addEventListener("click", () => openImageLightbox(imgEl.src));
        });
        galleryEl.querySelectorAll(".note-image-delete").forEach((btn) => {
          btn.addEventListener("click", async (e) => {
            e.stopPropagation();
            const imageId = btn.dataset.imageId;
            if (!await uiConfirm("이 이미지를 삭제할까요?")) return;
            try {
              const res = await fetch(`/api/sets/${encodeURIComponent(setId)}/verify/${entry.id}/images/${imageId}`, { method: "DELETE" });
              if (!res.ok) throw new Error(`HTTP ${res.status}`);
              entry.images = (entry.images || []).filter((img) => img.id !== imageId);
              refreshOnPage();
            } catch (err) {
              uiAlert(`삭제 실패: ${err.message}`);
            }
          });
        });
      }
    });
    container.querySelectorAll(".verify-delete-btn").forEach((btn) => {
      btn.addEventListener("click", async () => {
        const vid = btn.dataset.vid;
        if (!await uiConfirm("이 검증 기록을 삭제할까요?")) return;
        try {
          const res = await fetch(`/api/sets/${encodeURIComponent(setId)}/verify/${vid}`, { method: "DELETE" });
          if (!res.ok) throw new Error(`HTTP ${res.status}`);
          entries = entries.filter((e) => e.id !== vid);
          const s = state.sets.find((x) => x.set_id === setId);
          if (s) s.verifications = entries;
          refreshOnPage();
        } catch (err) {
          uiAlert(`삭제 실패: ${err.message}`);
        }
      });
    });
    container.querySelectorAll(".verify-add-image-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        container.querySelector(`.verify-image-file-input[data-vid="${btn.dataset.vid}"]`).click();
      });
    });
    container.querySelectorAll(".verify-image-file-input").forEach((input) => {
      input.addEventListener("change", async (e) => {
        const vid = input.dataset.vid;
        const entry = entries.find((x) => x.id === vid);
        for (const file of [...e.target.files]) {
          if (!file.type.startsWith("image/")) continue;
          try {
            const dataUrl = await new Promise((resolve, reject) => {
              const reader = new FileReader();
              reader.onload = () => resolve(reader.result);
              reader.onerror = reject;
              reader.readAsDataURL(file);
            });
            const base64 = dataUrl.split(",")[1];
            const res = await fetch(`/api/sets/${encodeURIComponent(setId)}/verify/${vid}/images`, {
              method: "POST", headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ content_type: file.type, data_base64: base64, caption: file.name }),
            });
            const data = await res.json();
            if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
            entry.images = entry.images || [];
            entry.images.push(data.image);
          } catch (err) {
            uiAlert(`이미지 업로드 실패: ${err.message}`);
          }
        }
        e.target.value = "";
        refreshOnPage();
      });
    });
  }

  refreshOnPage();

  currentVerifyTarget = {
    setId,
    setData,
    materialOptions,
    refreshOnPage,
    wireHistoryList,
    getEntries: () => entries,
  };

  $("#openVerifyModalBtn").addEventListener("click", () => {
    $("#vfMaterial").innerHTML = materialOptions();
    $("#vfDetected").value = "true";
    $("#vfCt").value = "";
    $("#vfNote").value = "";
    $("#vfFluorophore").value = setData.fluorophore || "";
    $("#vfSaveStatus").textContent = "";
    $("#vfInSilicoPreview").innerHTML = renderVerifyInSilicoPreview(setData, $("#vfMaterial").value);
    $("#vfHistoryList").innerHTML = renderVerificationHistory(entries, setId);
    wireHistoryList($("#vfHistoryList"));
    $("#verifyModal").hidden = false;
  });
}

function renderNoteImageGalleryFor(images, setId, vid, kind = "sets") {
  if (!images || !images.length) return "";
  return images.map((img) => `
    <div class="note-image-thumb" data-image-id="${img.id}">
      <img src="/api/${kind}/${encodeURIComponent(setId)}/verify/${vid}/images/${img.id}" alt="${img.caption || "검증 이미지"}" loading="lazy">
      <button class="note-image-delete" data-image-id="${img.id}" title="이미지 삭제">&times;</button>
      ${img.caption ? `<div class="note-image-caption">${img.caption}</div>` : ""}
    </div>`).join("");
}

$("#verifyModalClose").addEventListener("click", () => { $("#verifyModal").hidden = true; });
$("#verifyModal").addEventListener("click", (e) => { if (e.target.id === "verifyModal") $("#verifyModal").hidden = true; });

// #vfMaterial and #vfSaveBtn are singleton elements inside the shared #verifyModal (not
// recreated per render, unlike #openVerifyModalBtn which lives inside #detail's own template).
// Registering their listeners once here and reading currentVerifyTarget avoids stacking a new
// closure — bound to whichever set was open at the time — on every wireVerification() call.
$("#vfMaterial").addEventListener("change", () => {
  if (!currentVerifyTarget) return;
  $("#vfInSilicoPreview").innerHTML = renderVerifyInSilicoPreview(currentVerifyTarget.setData, $("#vfMaterial").value);
});

$("#vfSaveBtn").addEventListener("click", async () => {
  if (!currentVerifyTarget) return;
  const { setId, refreshOnPage, wireHistoryList, getEntries } = currentVerifyTarget;
  const statusEl = $("#vfSaveStatus");
  const material = $("#vfMaterial").value;
  if (!material) { statusEl.textContent = "표준물질을 선택하세요"; return; }
  statusEl.textContent = "저장 중...";
  try {
    const res = await fetch(`/api/sets/${encodeURIComponent(setId)}/verify`, {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({
        standard_material: material,
        detected: $("#vfDetected").value === "" ? "" : $("#vfDetected").value === "true",
        ct_value: $("#vfCt").value,
        note: $("#vfNote").value,
        fluorophore: $("#vfFluorophore").value,
      }),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    const entries = getEntries();
    entries.push(data.verification);
    const s = state.sets.find((x) => x.set_id === setId);
    if (s) s.verifications = entries;
    statusEl.textContent = "저장됨 — 아래에서 결과 이미지를 추가할 수 있습니다";
    $("#vfNote").value = "";
    $("#vfCt").value = "";
    $("#vfHistoryList").innerHTML = renderVerificationHistory(entries, setId);
    wireHistoryList($("#vfHistoryList"));
    refreshOnPage();
  } catch (err) {
    statusEl.textContent = `저장 실패: ${err.message}`;
  }
});

// ---------------- wet-lab verification for a GROUP (union of sets tested as one combo) ----------------
// Deliberately simpler than the per-set version above: a group has no single standard_genomes
// mismatch table to auto-derive an in-silico prediction from (member_details is the lightweight
// list-view shape, not the full per-genome detail), so this just records the wet-lab result
// itself (material, detected/not, Ct, photos) without the automatic in-silico-vs-wet-lab compare
// pill. Uses its own #grpVerifyModal (defined fresh inside renderGroupDetail's template each
// render) rather than the shared #verifyModal above, so its listeners don't pile up the way the
// per-set modal's do across repeated opens (that modal is a page-level singleton reused by every
// set; this one is destroyed and recreated with the rest of #detail every time).

function renderGroupVerificationHistory(entries, groupId) {
  if (!entries || !entries.length) return "<p class='note'>아직 등록된 실험 검증 기록이 없습니다.</p>";
  return entries.slice().reverse().map((entry) => {
    const wetLabel = entry.detected === true ? "검출(양성)" : entry.detected === false ? "미검출(음성)" : "판정보류";
    const wetCls = entry.detected === true ? "pill-ok" : entry.detected === false ? "pill-bad" : "";
    return `
    <div class="verify-entry" data-vid="${entry.id}">
      <div class="verify-entry-head">
        <span class="mono">${escHtml(entry.standard_material)}</span>
        <span class="pill ${wetCls}">${wetLabel}</span>
        ${entry.fluorophore ? `<span class="pill" title="형광물질(채널)">${escHtml(entry.fluorophore)}</span>` : ""}
        <span class="note">${escHtml(entry.created_at)}${entry.ct_value ? ` &middot; Ct ${escHtml(entry.ct_value)}` : ""}</span>
        <button class="btn btn-secondary grp-verify-delete-btn" data-vid="${entry.id}" title="검증 기록 삭제">삭제</button>
      </div>
      ${entry.note ? `<div class="note" style="margin:4px 0;">${escHtml(entry.note)}</div>` : ""}
      <div class="note-image-gallery" id="grpVhGallery-${entry.id}"></div>
      <div style="margin-top:6px;">
        <button class="btn btn-secondary grp-verify-add-image-btn" data-vid="${entry.id}">📷 결과 이미지 추가</button>
        <input type="file" class="grp-verify-image-file-input" data-vid="${entry.id}" accept="image/png,image/jpeg,image/gif,image/webp" multiple hidden>
      </div>
    </div>`;
  }).join("");
}

function wireGroupVerification(g) {
  const groupId = g.group_id;
  let entries = (g.verifications || []).slice();

  function refresh() {
    const target = $("#grpVerifyHistory");
    if (!target) return;
    target.innerHTML = renderGroupVerificationHistory(entries, groupId);
    wireHistoryList(target);
  }

  function wireHistoryList(container) {
    entries.forEach((entry) => {
      const galleryEl = container.querySelector(`#grpVhGallery-${entry.id}`);
      if (galleryEl) {
        galleryEl.innerHTML = renderNoteImageGalleryFor(entry.images, groupId, entry.id, "groups");
        galleryEl.querySelectorAll(".note-image-thumb img").forEach((imgEl) => {
          imgEl.addEventListener("click", () => openImageLightbox(imgEl.src));
        });
        galleryEl.querySelectorAll(".note-image-delete").forEach((btn) => {
          btn.addEventListener("click", async (e) => {
            e.stopPropagation();
            const imageId = btn.dataset.imageId;
            if (!await uiConfirm("이 이미지를 삭제할까요?")) return;
            try {
              const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}/verify/${entry.id}/images/${imageId}`, { method: "DELETE" });
              if (!res.ok) throw new Error(`HTTP ${res.status}`);
              entry.images = (entry.images || []).filter((img) => img.id !== imageId);
              refresh();
            } catch (err) {
              uiAlert(`삭제 실패: ${err.message}`);
            }
          });
        });
      }
    });
    container.querySelectorAll(".grp-verify-delete-btn").forEach((btn) => {
      btn.addEventListener("click", async () => {
        const vid = btn.dataset.vid;
        if (!await uiConfirm("이 검증 기록을 삭제할까요?")) return;
        try {
          const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}/verify/${vid}`, { method: "DELETE" });
          if (!res.ok) throw new Error(`HTTP ${res.status}`);
          entries = entries.filter((e) => e.id !== vid);
          refresh();
        } catch (err) {
          uiAlert(`삭제 실패: ${err.message}`);
        }
      });
    });
    container.querySelectorAll(".grp-verify-add-image-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        container.querySelector(`.grp-verify-image-file-input[data-vid="${btn.dataset.vid}"]`).click();
      });
    });
    container.querySelectorAll(".grp-verify-image-file-input").forEach((input) => {
      input.addEventListener("change", async (e) => {
        const vid = input.dataset.vid;
        const entry = entries.find((x) => x.id === vid);
        for (const file of [...e.target.files]) {
          if (!file.type.startsWith("image/")) continue;
          try {
            const dataUrl = await new Promise((resolve, reject) => {
              const reader = new FileReader();
              reader.onload = () => resolve(reader.result);
              reader.onerror = reject;
              reader.readAsDataURL(file);
            });
            const base64 = dataUrl.split(",")[1];
            const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}/verify/${vid}/images`, {
              method: "POST", headers: { "Content-Type": "application/json" },
              body: JSON.stringify({ content_type: file.type, data_base64: base64, caption: file.name }),
            });
            const data = await res.json();
            if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
            entry.images = entry.images || [];
            entry.images.push(data.image);
          } catch (err) {
            uiAlert(`이미지 업로드 실패: ${err.message}`);
          }
        }
        e.target.value = "";
        refresh();
      });
    });
  }

  refresh();

  $("#grpOpenVerifyBtn").addEventListener("click", () => {
    $("#grpVfDetected").value = "true";
    $("#grpVfCt").value = "";
    $("#grpVfNote").value = "";
    $("#grpVfFluorophore").value = "";
    $("#grpVfSaveStatus").textContent = "";
    $("#grpVerifyModal").hidden = false;
  });
  $("#grpVerifyModalClose").addEventListener("click", () => { $("#grpVerifyModal").hidden = true; });
  $("#grpVerifyModal").addEventListener("click", (e) => { if (e.target.id === "grpVerifyModal") $("#grpVerifyModal").hidden = true; });

  $("#grpVfSaveBtn").addEventListener("click", async () => {
    const statusEl = $("#grpVfSaveStatus");
    const material = $("#grpVfMaterial").value;
    if (!material) { statusEl.textContent = "표준물질을 선택하세요"; return; }
    statusEl.textContent = "저장 중...";
    try {
      const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}/verify`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({
          standard_material: material,
          detected: $("#grpVfDetected").value === "" ? "" : $("#grpVfDetected").value === "true",
          ct_value: $("#grpVfCt").value,
          note: $("#grpVfNote").value,
          fluorophore: $("#grpVfFluorophore").value,
        }),
      });
      const data = await res.json();
      if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
      entries.push(data.verification);
      statusEl.textContent = "저장됨";
      $("#grpVfNote").value = ""; $("#grpVfCt").value = "";
      refresh();
    } catch (err) {
      statusEl.textContent = `저장 실패: ${err.message}`;
    }
  });
}

// ---------------- phylogenetic tree (nextstrain auspice.json, rendered via PearTree) ----------------
// PearTree (artic-network/peartree, the FigTree successor) is a real embeddable tree viewer,
// self-hosted as vendor/peartree.bundle.min.js. That bundle injects ~4000 global CSS rules into
// whatever document loads it (it's built to own the whole page), so it must NEVER be loaded in
// the main document -- it's loaded only inside vendor/peartree_frame.html, mounted here as a
// same-origin <iframe>, which keeps its styling fully contained. Being same-origin lets us call
// straight into iframe.contentWindow.PearTreeEmbed and keep the controller object embed()
// returns (loadTree(), per PearTree's own embedded-api.md, lets us swap the tree without
// rebuilding the iframe -- used below to make the legend a real status filter).
//
// It reads Newick/NEXUS, not auspice v2 JSON directly, so we convert the fetched auspice tree
// to NEXUS here, encoding each tip's validation status/lineage/strain as a NEXUS metacomment.

state.auspiceCache = {};

async function loadAuspiceTree(seg) {
  const key = seg.toLowerCase();
  if (state.auspiceCache[key]) return state.auspiceCache[key];
  const data = await fetchJSON(`/api/auspice/${key}`);
  state.auspiceCache[key] = data;
  return data;
}

function statusOf(validatedIds, populationIds, gid) {
  if (validatedIds.has(gid)) return "validated";
  if (populationIds.has(gid)) return "notvalidated"; // in our curated DB, but this set/group didn't pass it
  return "absent"; // not in our curated DB at all (filtered out for length/N, etc.)
}

function collectTipStatusCounts(node, validatedIds, populationIds, counts) {
  if (!node.children || !node.children.length) {
    const accession = node.node_attrs?.accession || node.name;
    counts[statusOf(validatedIds, populationIds, accession)]++;
    counts.total++;
    return;
  }
  for (const c of node.children) collectTipStatusCounts(c, validatedIds, populationIds, counts);
}

function nexusSafeToken(text) {
  return String(text).replace(/[^A-Za-z0-9_.]/g, "_");
}

// PearTree's own categorical palette for a "status" trait picks its colours by itself, so the
// tree's colours never matched the legend buttons above it (green/red/grey). Its "user_colour"
// mode, however, uses each annotation value verbatim as the colour -- so every tip (and every
// internal branch whose tips all share one status) is annotated with exactly the legend's
// colour, and mixed internal branches get a neutral grey.
const TREE_MIXED_COLOR = "#7d8794";

function auspiceNodeToNewick(node, validatedIds, populationIds, parentDiv) {
  return _auspiceNode(node, validatedIds, populationIds, parentDiv).text;
}

function _auspiceNode(node, validatedIds, populationIds, parentDiv) {
  const div = typeof node.node_attrs?.div === "number" ? node.node_attrs.div : parentDiv;
  const bl = Math.max(0, div - parentDiv);
  const children = node.children || [];
  if (!children.length) {
    const accession = node.node_attrs?.accession || node.name;
    const status = statusOf(validatedIds, populationIds, accession);
    const lineage = node.node_attrs?.s_lineage?.value || node.node_attrs?.l_lineage?.value ||
      node.node_attrs?.gpc_lineage?.value || "NA";
    const strain = nexusSafeToken(node.name || accession);
    const tip = "T" + nexusSafeToken(accession);
    const color = TREE_STATUS_META[status].color;
    return { statuses: new Set([status]),
      text: `${tip}[&status=${status},lineage=${nexusSafeToken(lineage)},strain=${strain},user_colour="${color}"]:${bl.toFixed(6)}` };
  }
  const kids = children.map((c) => _auspiceNode(c, validatedIds, populationIds, div));
  const statuses = new Set();
  kids.forEach((k) => k.statuses.forEach((st) => statuses.add(st)));
  const color = statuses.size === 1 ? TREE_STATUS_META[[...statuses][0]].color : TREE_MIXED_COLOR;
  return { statuses, text: `(${kids.map((k) => k.text).join(",")})[&user_colour="${color}"]:${bl.toFixed(6)}` };
}

function auspiceToNexus(treeNode, validatedIds, populationIds) {
  const treeBody = auspiceNodeToNewick(treeNode, validatedIds, populationIds, 0);
  return `#NEXUS\nbegin trees;\n tree lasv_tree = [&R] ${treeBody};\nend;\n`;
}

function pruneAuspiceTree(node, keepFn) {
  // Drop tips keepFn rejects and collapse any internal node left with no children --
  // used to make the legend act as a real status filter on the rendered tree.
  if (!node.children || !node.children.length) {
    const accession = node.node_attrs?.accession || node.name;
    return keepFn(accession) ? node : null;
  }
  const kids = node.children.map((c) => pruneAuspiceTree(c, keepFn)).filter(Boolean);
  if (!kids.length) return null;
  return { ...node, children: kids };
}

const TREE_STATUS_META = {
  validated: { color: "#54A24B", label: "검증 통과" },
  notvalidated: { color: "#E45756", label: "검증 실패/불일치" },
  absent: { color: "#9aa3af", label: "population DB에 없음" },
};

async function renderTreeSection(container, data, seg) {
  container.innerHTML = "<p class='note'>Nextstrain 계통수 불러오는 중... (처음 한 번은 다운로드 때문에 몇 초 걸릴 수 있습니다)</p>";
  try {
    const auspiceData = await loadAuspiceTree(seg);
    const validatedIds = new Set(data.genomes
      ? data.genomes.filter(isValidatedGenome).map((g) => g.genome_id)
      : (data.validated_ids || []));
    const populationIds = new Set(data.population_ids || []);
    const counts = { validated: 0, notvalidated: 0, absent: 0, total: 0 };
    collectTipStatusCounts(auspiceData.tree, validatedIds, populationIds, counts);

    const activeStatuses = new Set(Object.keys(TREE_STATUS_META));
    const legendHtml = Object.entries(TREE_STATUS_META).map(([key, m]) => `
      <button type="button" class="tree-legend-btn active" data-status="${key}">
        <span class="dot" style="background:${m.color}"></span>${m.label} (${counts[key]})
      </button>`).join("");

    container.innerHTML = `
      <div class="tree-legend legend" style="margin-bottom:8px;">
        ${legendHtml}
        <span>Nextstrain 원본 트리 tip 전체: ${counts.total}개</span>
        <span class="note">위 항목을 클릭하면 해당 상태의 tip만 트리에 남습니다 (다시 클릭하면 복귀)</span>
      </div>
      <div class="tree-outer peartree-outer"></div>
    `;
    const mount = container.querySelector(".tree-outer");
    const iframe = document.createElement("iframe");
    iframe.className = "peartree-frame";
    iframe.style.cssText = "width:100%;height:820px;border:none;display:block;";
    iframe.src = "vendor/peartree_frame.html";
    mount.appendChild(iframe);

    // Colour by the per-node "user_colour" annotation written in auspiceNodeToNewick() -- the
    // exact legend colours. (The settings object uses camelCase keys, not the dash-case names of
    // PearTree's own form controls.)
    const STATUS_COLOR_SETTINGS = {
      tipColourBy: "user_colour",
      branchColourBy: "user_colour",
      nodeColourBy: "user_colour",
      labelColourBy: "user_colour",
    };

    // embed() (and the controller's own loadTree) resolve/act async, so guard against a
    // legend click landing while a previous load is still in flight -- otherwise embed()
    // gets called twice and PearTree mounts a second instance inside the same container.
    let controller = null;
    let busy = false;
    let pending = false;
    const loadFiltered = async () => {
      if (busy) { pending = true; return; }
      busy = true;
      try {
        const pruned = pruneAuspiceTree(auspiceData.tree, (acc) => activeStatuses.has(statusOf(validatedIds, populationIds, acc)));
        const nexus = pruned ? auspiceToNexus(pruned, validatedIds, populationIds) : "#NEXUS\nbegin trees;\n tree empty = [&R] ();\nend;\n";
        if (controller) {
          await controller.loadTree(nexus, "tree.nex");
          await controller.applySettings?.(STATUS_COLOR_SETTINGS);
        } else if (iframe.contentWindow?.PearTreeEmbed) {
          controller = await iframe.contentWindow.PearTreeEmbed.embed({
            container: "mount", tree: nexus, height: "100%", settings: STATUS_COLOR_SETTINGS,
          });
          await controller.applySettings?.(STATUS_COLOR_SETTINGS);
        }
      } finally {
        busy = false;
        if (pending) { pending = false; loadFiltered(); }
      }
    };

    iframe.addEventListener("load", () => {
      if (!iframe.contentWindow?.PearTreeEmbed) {
        mount.innerHTML = "<p class='note'>PearTree 라이브러리를 불러오지 못했습니다.</p>";
        return;
      }
      loadFiltered().catch((e) => { mount.innerHTML = `<p class='note'>PearTree 로드 실패: ${e.message}</p>`; });
    }, { once: true });

    container.querySelectorAll(".tree-legend-btn").forEach((btn) => {
      btn.addEventListener("click", () => {
        const status = btn.dataset.status;
        if (activeStatuses.has(status)) {
          if (activeStatuses.size === 1) return; // keep at least one status visible
          activeStatuses.delete(status);
          btn.classList.remove("active");
        } else {
          activeStatuses.add(status);
          btn.classList.add("active");
        }
        loadFiltered();
      });
    });
  } catch (e) {
    container.innerHTML = `<p class="note">계통수 불러오기 실패: ${e.message}</p>`;
  }
}

// Where a pop-up (dialog, detail window, lightbox) is attached: while a card is shown fullscreen only that
// card's subtree is displayed, so a pop-up appended to <body> would open invisibly behind it.
function modalHost() {
  return document.fullscreenElement || document.body;
}

function toggleFullscreen(el) {
  if (document.fullscreenElement) {
    document.exitFullscreen();
  } else {
    (el.requestFullscreen || el.webkitRequestFullscreen)?.call(el);
  }
}

// ---------------------------------------------------------------------------------------
// Minimize/maximize controls for every result panel (.card and its chart/stat-card
// siblings) -- injected generically rather than hand-wired into each of the dozens of
// places that build one of these panels. A MutationObserver on the content area means any
// future card, from any view, gets the same controls automatically with zero extra code at
// its own call site.
const CARD_CONTROL_SELECTOR = ".card, .diversity-chart-card, .meta-map-card, .meta-stat-card";

// Collapsed state persists across reloads, keyed by the card's own title text -- a global
// "I don't need to see the X card by default" preference rather than a per-instance one, so
// collapsing e.g. "Population 커버리지" once applies to every set's detail page, not just
// the one open right now (matches how most apps remember a collapsed section as a standing
// preference, not transient scroll state).
const CARD_STATE_KEY = "lasv_card_collapsed_v1";
function loadCardCollapsedState() {
  try { return JSON.parse(localStorage.getItem(CARD_STATE_KEY) || "{}"); } catch { return {}; }
}
function saveCardCollapsedState(state) {
  try { localStorage.setItem(CARD_STATE_KEY, JSON.stringify(state)); } catch { /* private mode etc. -- just skip persisting */ }
}
function cardTitleKey(card, ownHeader) {
  const h3 = ownHeader?.querySelector("h3") || card.querySelector(":scope > h3");
  return h3 ? h3.textContent.trim() : null;
}

// How tall a collapsed card should stay (enough to show just its own header) -- measured
// per-card instead of one fixed pixel value, since a card's title line can wrap or run a
// different font-size than another's. Falls back to a fixed default for header-less cards,
// where the "first child" would otherwise BE the full content block.
function peekHeight(card, headerEl) {
  const ref = headerEl || card.querySelector(":scope > h3");
  if (!ref) return 48;
  const pad = parseFloat(getComputedStyle(card).paddingTop) || 16;
  const h = Math.ceil(ref.getBoundingClientRect().height + pad * 2);
  return h > 0 && h < 120 ? h : 48;
}

function setCardCollapsed(card, collapsed, headerEl) {
  card.classList.toggle("card-collapsed", collapsed);
  if (collapsed) card.style.setProperty("--card-peek-h", peekHeight(card, headerEl) + "px");
  const key = cardTitleKey(card, headerEl);
  if (key) {
    const state = loadCardCollapsedState();
    if (collapsed) state[key] = true; else delete state[key];
    saveCardCollapsedState(state);
  }
}

function injectCardControls(card) {
  if (card.dataset.ctrlsInjected) return;
  card.dataset.ctrlsInjected = "1";
  const savedState = loadCardCollapsedState();

  // .msa-card / .tree-card already have their own header row with a fullscreen button --
  // add just a matching minimize button into that row instead of a second, floating one.
  const ownHeader = card.querySelector(":scope > .msa-card-header");
  const ownFsBtn = ownHeader?.querySelector(".btn-fullscreen");
  if (ownHeader && ownFsBtn) {
    const minBtn = document.createElement("button");
    minBtn.className = "btn btn-secondary btn-fullscreen";
    minBtn.title = "최소화 / 펼치기";
    minBtn.innerHTML = "&#8722;";
    minBtn.addEventListener("click", (e) => {
      e.stopPropagation();
      const collapsed = !card.classList.contains("card-collapsed");
      setCardCollapsed(card, collapsed, ownHeader);
      minBtn.innerHTML = collapsed ? "&#43;" : "&#8722;";
    });
    ownHeader.insertBefore(minBtn, ownFsBtn);
    if (savedState[cardTitleKey(card, ownHeader)]) { setCardCollapsed(card, true, ownHeader); minBtn.innerHTML = "&#43;"; }
    return;
  }

  const controls = document.createElement("div");
  controls.className = "card-controls";
  const minBtn = document.createElement("button");
  minBtn.className = "card-ctrl-btn"; minBtn.title = "최소화 / 펼치기"; minBtn.innerHTML = "&#8722;";
  minBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    const collapsed = !card.classList.contains("card-collapsed");
    setCardCollapsed(card, collapsed);
    minBtn.innerHTML = collapsed ? "&#43;" : "&#8722;";
  });
  const maxBtn = document.createElement("button");
  maxBtn.className = "card-ctrl-btn"; maxBtn.title = "최대화 (전체화면)"; maxBtn.innerHTML = "&#9974;";
  maxBtn.addEventListener("click", (e) => {
    e.stopPropagation();
    setCardCollapsed(card, false);
    toggleFullscreen(card);
  });
  controls.appendChild(minBtn);
  controls.appendChild(maxBtn);
  card.appendChild(controls);
  if (savedState[cardTitleKey(card, null)]) { setCardCollapsed(card, true); minBtn.innerHTML = "&#43;"; }
}

function attachCardControls(root) {
  root.querySelectorAll(CARD_CONTROL_SELECTOR).forEach(injectCardControls);
}

// attachCardControls()/translateNode() both mutate the very subtree this observer watches
// (appending control buttons, replacing translated text) -- without disconnecting first,
// each of those mutations would re-trigger this same callback, which mutates again, forever.
// disconnect() before, observe() again after is the standard-safe pattern for an observer
// whose own callback writes to its observed tree.
const cardControlsObserver = new MutationObserver((mutations) => {
  for (const m of mutations) {
    if (m.addedNodes.length) {
      const main = document.querySelector(".main");
      cardControlsObserver.disconnect();
      attachCardControls(main);
      translateNode(main);
      cardControlsObserver.observe(main, { childList: true, subtree: true });
      return;
    }
  }
});
cardControlsObserver.observe(document.querySelector(".main"), { childList: true, subtree: true });

// SeaLion (https://artic-network.github.io/sealion) is a standalone large-alignment
// viewer, not an embeddable library -- so instead of reimplementing it, export a FASTA
// the user can drag straight into it for a deeper look.
function exportMsaFasta(setData) {
  const msa = setData.msa;
  if (!msa || !msa.columns) return;
  const undot = (dots) => dots.split("").map((c, i) => (c === "." ? msa.consensus[i] : c)).join("");
  const lines = [];
  for (const key of ["forward", "probe", "reverse"]) {
    const rr = msa.ref_rows && msa.ref_rows[key];
    if (rr) lines.push(`>REFERENCE_${key.toUpperCase()}_${rr.strand === "-" ? "revcomp" : "fwd"}\n${rr.aligned}`);
  }
  lines.push(`>CONSENSUS\n${msa.consensus}`);
  for (const r of msa.rows) {
    lines.push(`>${r.genome_id}|${r.lineage}|${r.country || "?"}\n${undot(r.dots)}`);
  }
  const blob = new Blob([lines.join("\n") + "\n"], { type: "text/plain" });
  const url = URL.createObjectURL(blob);
  const a = document.createElement("a");
  a.href = url;
  a.download = `${setData.set_id}_alignment.fasta`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(url);
  window.open("https://artic-network.github.io/sealion", "_blank");
}

async function selectSet(setId) {
  state.activeId = setId;
  state.activeType = "set";
  state.msaSearch = "";
  state.msaFilterMode = "all";
  state.msaColFilters = {};
  renderSidebar();
  showMainView("detail");
  const detail = $("#detail");
  detail.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    const data = await fetchJSON(`/api/sets/${encodeURIComponent(setId)}`);
    renderDetail(data);
    if (data.host_scope) detail.insertAdjacentHTML("afterbegin", hostScopeBadgeHtml({ [(data.segment || "").toUpperCase()]: data.host_scope }));
  } catch (e) {
    detail.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
  }
}

async function selectGroup(groupId) {
  state.activeId = groupId;
  state.activeType = "group";
  renderSidebar();
  showMainView("detail");
  const detail = $("#detail");
  detail.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    const data = await fetchJSON(`/api/groups/${encodeURIComponent(groupId)}`);
    renderGroupDetail(data);
  } catch (e) {
    detail.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
  }
}

function renderSizeFilterNote(sizeFilter) {
  if (!sizeFilter) return "";
  const removed = (sizeFilter.combos_before_size_filter || 0) - (sizeFilter.combos_after_size_filter || 0);
  const anchorText = sizeFilter.anchor_size !== null && sizeFilter.anchor_size !== undefined
    ? `${Math.round(sizeFilter.anchor_size)}bp` : "정보 없음";
  const perSet = sizeFilter.per_set_sizes || {};
  const sourceLabel = { refgenome: "참조유전체 기준(정확)", population_median: "population 중앙값(근사)" };
  const perSetHtml = Object.keys(perSet).length ? `
    <div class="size-filter-detail">
      ${Object.entries(perSet).map(([sid, info]) => `<span class="mono">${sid}: ${Math.round(info.size)}bp (${sourceLabel[info.source] || info.source})</span>`).join("")}
    </div>` : "";
  return `
    <div class="note size-filter-note" style="margin-bottom:10px;">
      📏 amplicon 크기 필터 적용됨: 각 세트가 자기 자신과 정상적으로 결합할 때의 실제 amplicon 크기를
      (가능하면 참조유전체로 검증된 정확한 크기를, 없는 세트만 population 중앙값으로) 먼저 확인한 뒤,
      그중 가장 큰 값(${anchorText})의 2배(최대 ${sizeFilter.ideal_max_size}bp)를 넘는 조합은 실제
      PCR에서 증폭될 가능성이 낮다고 보고 제외했습니다${removed > 0 ? ` (${removed}개 조합 제외됨)` : ""}.
      두 primer가 기하학적으로는 짝을 이룰 수 있어도 너무 멀리 떨어져 있으면 실제 위협이 아니므로, 이
      그룹 자체의 정상 결합 크기를 기준으로 판단합니다.
      ${perSetHtml}
    </div>`;
}

function renderComboSummary(comboGenomes) {
  const byCombo = {};
  for (const r of comboGenomes) {
    const key = r.combo_label || r.combo_id;
    if (!byCombo[key]) byCombo[key] = { count: 0, pass: 0, isSameSet: r.is_same_set };
    byCombo[key].count++;
    if (r.p_mismatch != null && r.p_mismatch !== "" && verdictFor(r).ok) byCombo[key].pass++;      // F/R <= 1, 3' ends intact, probe perfect
  }
  const combos = Object.keys(byCombo).sort((a, b) => a.localeCompare(b));
  if (!combos.length) return "<p class='note'>발견된 결합 조합이 없습니다.</p>";
  return combos.map((c) => {
    const info = byCombo[c];
    const cls = info.isSameSet ? "combo-chip-ok" : "combo-chip-warn";
    const tag = info.isSameSet ? "정상 (같은 세트)" : "⚠ 교차 결합 (다른 세트)";
    return `<div class="combo-chip ${cls}"><span class="combo-chip-name mono">${c}</span><span class="combo-chip-meta">${tag} · ${info.count}개 게놈 결합 · <b>${info.pass}개 통과</b></span></div>`;
  }).join("");
}

function comboRowStatus(r) {
  // One combined, at-a-glance verdict per row -- mismatch-based pass/fail (verdictFor)
  // plus whether the binding itself looks thermodynamically/size-wise plausible
  // (assessBindingQuality) -- so a row can look "no significant mismatch" and still be
  // flagged if e.g. its Tm/ΔG is too weak to realistically amplify.
  const verdict = verdictFor(r);
  if (!verdict.ok) return nearBinding(r) ? { level: "warn", icon: "⚠️", label: "허용 초과" } : { level: "bad", icon: "🚫", label: "문제" };
  const quality = assessBindingQuality(r);
  if (quality.some((q) => q.ok === false)) return { level: "warn", icon: "⚠️", label: "주의" };
  return { level: "ok", icon: "✅", label: "정상" };
}

function renderComboGenomeTable(comboGenomes, setsById, idPrefix) {
  idPrefix = idPrefix || "combo";
  const search = ($(`#${idPrefix}Search`)?.value || "").toLowerCase();
  const crossOnly = state[`${idPrefix}CrossOnly`];
  const problemsOnly = state[`${idPrefix}ProblemsOnly`];
  let rows = comboGenomes.filter((r) => {
    if (crossOnly && r.is_same_set) return false;
    if (problemsOnly && comboRowStatus(r).level === "ok") return false;
    if (search) {
      const hay = `${r.genome_id} ${r.combo_id} ${r.combo_label || ""} ${r.lineage} ${r.country}`.toLowerCase();
      if (!hay.includes(search)) return false;
    }
    return true;
  });
  rows = rows.slice().sort((a, b) => a.genome_id.localeCompare(b.genome_id) || a.combo_id.localeCompare(b.combo_id));
  const matchedCount = rows.length;
  const COMBO_ROW_CAP = 500;
  const truncated = rows.length > COMBO_ROW_CAP;
  if (truncated) rows = rows.slice(0, COMBO_ROW_CAP);

  const tbody = $(`#${idPrefix}TableBody`);
  if (!tbody) return;
  tbody.innerHTML = "";
  const frag = document.createDocumentFragment();
  for (const r of rows) {
    const status = comboRowStatus(r);
    const tr = el("tr", "genome-row" + (r.is_same_set ? "" : " problem-row"));
    tr.dataset.genomeId = r.genome_id;
    tr.innerHTML = `
      <td><span class="combo-status combo-status-${status.level}">${status.icon} ${status.label}</span></td>
      <td class="mono">${r.genome_id}</td>
      <td>${r.lineage || "NA"}</td>
      <td class="mono combo-id-cell">${r.combo_label || r.combo_id}${r.is_same_set ? "" : " <span class='combo-warn-tag'>⚠ 교차</span>"}</td>
      <td><span class="msa-mm-tag ${mmClass(r.f_mismatch)}">${mmText(r.f_mismatch)}</span></td>
      <td>${r.p_set_id ? `<span class="msa-mm-tag ${mmClass(r.p_mismatch)}">${mmText(r.p_mismatch)}</span>` : "<span class='note'>—</span>"}</td>
      <td><span class="msa-mm-tag ${mmClass(r.r_mismatch)}">${mmText(r.r_mismatch)}</span></td>
      <td class="mono">${r.amplicon_start}-${r.amplicon_end}</td>`;
    // f_oligo_id/r_oligo_id are "<set_id>::F" or "<set_id>::R" -- the oligo that ends up
    // playing the geometric left(fp)/right(rp) role isn't always that set's own F/R
    // (see group_cross_binding.py), so which sequence to actually display (and how to
    // label it) must be read from the resolved role, not assumed from f_set_id/r_set_id
    // alone -- using the wrong one here previously showed a completely different
    // sequence than what the combo_id / table row said was matched.
    const fRole = (r.f_oligo_id || "").split(ID_SEP).pop() || "F";
    const rRole = (r.r_oligo_id || "").split(ID_SEP).pop() || "R";
    const seqForRole = (setId, role) => {
      const s = setsById[setId];
      if (!s) return "";
      return role === "R" ? (s.R_seq || "") : (s.F_seq || "");
    };
    const setData = {
      F_seq: seqForRole(r.f_set_id, fRole),
      R_seq: seqForRole(r.r_set_id, rRole),
      P_seq: r.p_set_id ? (setsById[r.p_set_id]?.P_seq || "") : "",
    };
    const labels = {
      f: `${r.f_set_name || r.f_set_id} · ${fRole}`,
      p: r.p_set_id ? `${r.p_set_name || r.p_set_id} · P` : "PROBE",
      r: `${r.r_set_name || r.r_set_id} · ${rRole}`,
    };
    tr.addEventListener("click", () => openGenomeModal(r, setData, labels, r.combo_label || r.combo_id));
    frag.appendChild(tr);
  }
  tbody.appendChild(frag);
  const countLabel = $(`#${idPrefix}TableCount`);
  if (countLabel) {
    countLabel.textContent = truncated
      ? `${rows.length} / ${matchedCount}개 표시 (검색으로 좁혀보세요, 전체 ${comboGenomes.length}행)`
      : `${matchedCount} / ${comboGenomes.length}개 표시`;
  }
}

function renderComboSection(idPrefix, title, desc, extraHtml) {
  return `
    <div class="card" style="margin-bottom:20px;">
      <h3>${title}</h3>
      <div class="note" style="margin-bottom:10px;">${desc}</div>
      ${extraHtml || ""}
      <div id="${idPrefix}Summary" class="combo-summary"></div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <div class="msa-card-header">
        <h3 style="margin:0;">게놈별 ${title} F / Probe / R 미스매치</h3>
        <span id="${idPrefix}TableCount" class="note"></span>
      </div>
      <div class="table-toolbar">
        <input id="${idPrefix}Search" type="text" placeholder="ID, 조합, 리니지 검색...">
        <label><input id="${idPrefix}CrossOnlyChk" type="checkbox"> 교차 결합만 보기</label>
        <label><input id="${idPrefix}ProblemsOnlyChk" type="checkbox"> 문제/주의만 보기</label>
      </div>
      <div class="genome-table-wrap">
        <table class="genome-table">
          <thead>
            <tr>
              <th>상태</th><th>ID</th><th>Lineage</th><th>조합 (F set x R set)</th>
              <th>F mismatch</th><th>Probe mismatch</th><th>R mismatch</th><th>Amplicon</th>
            </tr>
          </thead>
          <tbody id="${idPrefix}TableBody"></tbody>
        </table>
      </div>
      <div class="legend" style="margin-top:8px;">
        <span><span class="combo-status combo-status-ok">✅ 정상</span></span>
        <span><span class="combo-status combo-status-warn">⚠️ 주의</span></span>
        <span><span class="combo-status combo-status-bad">🚫 문제</span></span>
      </div>
    </div>`;
}

function wireComboSection(idPrefix, comboGenomes, setsById) {
  state[`${idPrefix}CrossOnly`] = false;
  state[`${idPrefix}ProblemsOnly`] = false;
  $(`#${idPrefix}Summary`).innerHTML = renderComboSummary(comboGenomes);
  renderComboGenomeTable(comboGenomes, setsById, idPrefix);
  $(`#${idPrefix}Search`).addEventListener("input", () => renderComboGenomeTable(comboGenomes, setsById, idPrefix));
  $(`#${idPrefix}CrossOnlyChk`).addEventListener("change", (e) => {
    state[`${idPrefix}CrossOnly`] = e.target.checked;
    renderComboGenomeTable(comboGenomes, setsById, idPrefix);
  });
  $(`#${idPrefix}ProblemsOnlyChk`).addEventListener("change", (e) => {
    state[`${idPrefix}ProblemsOnly`] = e.target.checked;
    renderComboGenomeTable(comboGenomes, setsById, idPrefix);
  });
}

function renderGroupMemberCoverage(memberDetails) {
  if (!memberDetails || !memberDetails.length) return "<p class='note'>세트 정보가 없습니다.</p>";
  const rows = memberDetails.map((s) => {
    const sm = s.summary || {};
    return `<tr>
      <td class="mono">${escHtml(s.display_name || s.set_id)}</td>
      <td>${fmtPct(sm.primer_coverage_mm0_pct)}</td>
      <td>${fmtPct(sm.primer_coverage_mm1_pct)}</td>
      <td>${fmtPct(sm.primer_coverage_mm2_pct)}</td>
      <td>${fmtPct(sm.probe_coverage_mm0_pct)}</td>
      <td style="font-weight:700;color:var(--accent);">${fmtPct(sm.validated_pct)}</td>
    </tr>`;
  }).join("");
  return `
    <div class="genome-table-wrap" style="max-height:none;">
      <table class="genome-table">
        <thead><tr>
          <th>프라이머 세트</th><th>primer mm=0</th><th>primer mm=1</th><th>primer mm=2</th>
          <th>probe mm=0</th><th>validated</th>
        </tr></thead>
        <tbody>${rows}</tbody>
      </table>
    </div>`;
}

// Feeds the exact same narrateSetAssessment()/build_assessment() heuristic a group's
// union-level data -- per-set-only checks (probe dropout, single amplicon size,
// single-genome binding) are left out via empty/absent fields (both narrateSetAssessment
// and its Python twin already skip a check gracefully when its inputs are missing) rather
// than duplicating the scoring logic with a parallel group-specific heuristic. Mirrors
// webapp/server.py's build_group_evaluation() so the badge shown here matches what
// `g.evaluation` (computed server-side, e.g. for a future list view) would say.
function groupAssessmentInput(g) {
  const members = g.member_details || [];
  const stdOk = members.some((m) => (m.standard || {}).ok);
  const matched = Array.from(new Set(members.filter((m) => (m.standard || {}).ok)
    .flatMap((m) => m.standard.matched || [])));
  return {
    // coverage of the group = union of the sets' own results plus the genomes passed by a cross-binding pairing
    summary: { validated_pct: g.summary.total_genomes && g.validated_ids_with_cross
      ? (g.validated_ids_with_cross.length / g.summary.total_genomes) * 100 : g.summary.validated_pct },
    standard: { ok: stdOk, matched },
    gene: {},
    binding: [],
    dimer_hairpin: g.dimer_hairpin,
    cross_binding: g.combo_genomes || [],
  };
}

function renderGroupDetail(g) {
  const detail = $("#detail");
  state.detailGroup = g;          // read by the group analyze / re-analyze buttons
  const totals = g.lineage_totals || {};
  // every lineage of the population (also those the group covers with 0 genomes -- a gap worth seeing)
  const lineageSource = g.lineage_counts_with_cross || g.lineage_counts;     // union of the sets' own results + cross-binding passes
  const lineageCount = (l) => lineageSource[l] || 0;
  const lineages = Array.from(new Set([...Object.keys(totals), ...Object.keys(lineageSource)])).sort((a, b) => lineageCount(b) - lineageCount(a));
  const lineageHtml = lineages.map((l) => {
    const denom = totals[l] || 0;
    const pct = denom ? ((lineageCount(l) / denom) * 100).toFixed(0) : null;
    return `
    <div class="cov-row" style="grid-template-columns:110px 1fr 110px;cursor:default;">
      <div class="label">${l}</div>
      <div class="cov-bar-track"><div class="cov-bar-fill validated" style="width:${denom ? Math.min(100, (lineageCount(l) / denom) * 100) : 0}%"></div></div>
      <div class="cov-val" style="margin-left:0;">${lineageCount(l)}${denom ? ` / ${denom} (${pct}%)` : ""}</div>
    </div>`;
  }).join("");
  const memberChips = g.members.map((m) => `<span class="group-chip" data-set-id="${m}">${escHtml((g.member_names || {})[m] || m)}</span>`).join("");

  detail.innerHTML = `
    <div class="detail-header">
      <h2>&#9633; ${g.label}</h2>
      <button id="renameGroupBtn" class="btn btn-secondary" title="그룹 이름 변경">&#9998; 이름 변경</button>
      <button id="analyzeGroupBtn" class="btn btn-secondary" title="이 그룹만 분석합니다 (합집합·교차 결합·표준물질/참조유전체 조합). 이미 최신인 단계는 건너뜁니다.">&#9654; 이 그룹 분석</button>
      <button id="reanalyzeGroupBtn" class="btn btn-secondary" title="입력이 같아도 이 그룹의 조합 결과를 처음부터 다시 계산합니다 (구성 세트의 개별 결과는 건드리지 않음).">&#8635; 이 그룹 재분석</button>
      <span class="pill pill-${g.segment}">${g.segment}</span>
      <span class="pill">ALL (${g.members.length}개 세트 합집합)</span>
      <button id="deleteGroupBtn" class="btn btn-secondary" style="margin-left:auto;color:var(--bad);border-color:var(--bad);">그룹 삭제</button>
    </div>
    ${finalResultBoxHtml(g.group_id)}
    <div class="detail-sub">구성 세트: ${memberChips}</div>

    ${summaryWithCurve(renderSetAssessment(groupAssessmentInput(g)), g.segment, g.members)}

    <div class="card-grid">
      <div class="card">
        <h3>그룹 커버리지</h3>
        <div class="metric-list">
          <div class="row" style="padding-top:2px;">
            <span><b>합집합 (세트별 단독)</b></span>
            <span class="v" style="color:var(--accent);">${g.summary.validated} / ${g.summary.total_genomes} (${fmtPct(g.summary.validated_pct)})</span>
          </div>
          <div class="row">
            <span><b>교차 결합 포함</b></span>
            <span class="v" style="color:var(--accent);">${g.summary.validated_with_cross ?? g.summary.validated} / ${g.summary.total_genomes} (${fmtPct(g.summary.validated_with_cross_pct ?? g.summary.validated_pct)})</span>
          </div>
        </div>
        <div class="note" style="margin-top:6px;">교차 결합 포함: 세트 간 조합(예: 세트2의 F/R + 세트1의 프로브)이 F/R mismatch&le;1(3′ 끝 일치), 프로브 완전일치로 통과하는 게놈까지 센 수입니다. 개별 세트에서는 안 나와도 그룹으로 묶으면 통과하는 게놈이 여기 나타납니다.</div>
        <div class="note" style="margin-top:8px;">이 그룹의 세트 중 하나라도 통과(validated)하면 그 게놈은 합집합에 포함됩니다.</div>
      </div>
      <div class="card">
        <h3>리니지 분포</h3>
        <div class="note" style="margin-bottom:6px;">교차 결합 포함 통과 게놈 기준 (${(g.validated_ids_with_cross || g.validated_ids || []).length}개)</div>
        <div class="cov-rows" style="margin-bottom:0;">${lineageHtml || "<p class='note'>없음</p>"}</div>
      </div>
      <div class="card">
        <h3>이차구조 (Self/Cross-dimer &middot; Hairpin, 그룹 전체)</h3>
        <div class="note" style="margin-bottom:6px;">구성 세트 전체의 F/Probe/R을 모두 모아 검사 — 한 세트의 올리고가 <b>다른 세트</b>의 올리고와 결합하는 경우까지 포함합니다.</div>
        ${renderDimerHairpin(g.dimer_hairpin)}
      </div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>세트별 Population 커버리지</h3>
      ${renderGroupMemberCoverage(g.member_details)}
    </div>

    ${renderSetDiversitySection(g)}

    <div class="card tree-card" style="margin-bottom:20px;">
      <div class="msa-card-header">
        <h3 style="margin:0;">Nextstrain 계통수 대비 커버리지</h3>
        <button id="treeFullscreenBtn" class="btn btn-secondary btn-fullscreen" title="전체화면">&#9974; 전체화면</button>
      </div>
      <div id="treeContainer">
        <button id="loadTreeBtn" class="btn btn-secondary">계통수 불러오기</button>
      </div>
    </div>

    ${renderComboSection("combo", "교차 결합 분석 (Cross-binding, Population)",
      `각 세트를 독립적으로 분석해 단순 합집합만 내는 대신, 이 그룹의 모든 F/R/Probe를 하나로 모아
      다시 mfeprimer 검색을 수행해 <b>세트 간에 의도치 않게 결합하는 조합</b>(예: 세트1의 F가 세트2의
      R과 결합, 세트2의 Probe가 세트1의 amplicon 안에 위치)까지 실제로 찾아낸 결과입니다.`,
      renderSizeFilterNote(g.size_filter))}

    ${renderComboSection("standardCombo", "표준물질(KDCA) 매치 (Cross-binding)",
      `population 대신 KDCA 표준물질 DB에 대해 같은 방식(그룹의 모든 F/R/Probe를 한 번에 검색)으로
      교차 결합까지 포함해 확인합니다 — 세트 하나만 봤을 때는 안 보이던, 다른 세트의 Probe가 표준물질에
      결합하는 경우까지 여기서 드러납니다.`)}

    ${renderComboSection("refgenomeCombo", "참조유전체 결합위치 (Cross-binding)",
      `단일 참조유전체(S=AY628203, L=AY628202)에 대해 같은 방식으로 교차 결합까지 포함해 결합위치를
      확인합니다.`)}

    <div class="card" style="margin-bottom:20px;">
      <h3 style="display:flex;align-items:center;justify-content:space-between;gap:6px;">
        🧪 실험(Wet-lab) 검증 (그룹 전체)
        <button id="grpOpenVerifyBtn" class="btn btn-primary">+ 검증 결과 추가</button>
      </h3>
      <div class="note" style="margin-bottom:10px;">
        이 그룹(합집합)을 하나의 조합으로 실험한 결과를 기록합니다. 개별 세트만 테스트한 결과는 해당
        세트 상세 페이지에 기록하세요.
      </div>
      <div id="grpVerifyHistory"></div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>메모 / 노트 (그룹 전체)</h3>
      <textarea id="setNoteInput" rows="3" placeholder="이 그룹(조합)에 대한 메모를 남겨두세요... 사진은 여기에 붙여넣기(Ctrl/Cmd+V)하거나 아래 버튼으로 추가할 수 있습니다." style="width:100%;resize:vertical;padding:8px 10px;border-radius:8px;border:1px solid var(--border);background:var(--bg);color:var(--fg);font-size:12.5px;">${g.note || ""}</textarea>
      <div style="display:flex;align-items:center;gap:10px;margin-top:6px;">
        <button id="saveNoteBtn" class="btn btn-primary">저장</button>
        <button id="addNoteImageBtn" class="btn btn-secondary">📷 사진 추가</button>
        <input id="noteImageFileInput" type="file" accept="image/png,image/jpeg,image/gif,image/webp" multiple hidden>
        <span id="noteSaveStatus" class="note"></span>
      </div>
      <div id="noteImageGallery" class="note-image-gallery"></div>
    </div>

    <div id="grpVerifyModal" class="modal-overlay" hidden>
      <div class="modal-box verify-modal-box">
        <div class="modal-head">
          <div>
            <div class="modal-title">실험(Wet-lab) 검증 결과 추가 — 그룹 전체</div>
            <div class="modal-sub">표준물질(양성대조물질)로 이 그룹(구성 세트 전체)을 테스트한 결과를 기록합니다</div>
          </div>
          <button id="grpVerifyModalClose" class="btn btn-secondary">닫기</button>
        </div>
        <div class="verify-modal-body">
          <div class="form-grid">
            <label>표준물질 선택
              <select id="grpVfMaterial">${(g.standard_materials || []).map((m) => `<option value="${escHtml(m)}">${escHtml(m)}</option>`).join("")}</select>
            </label>
            <label>검출 여부 (실험 결과)
              <select id="grpVfDetected">
                <option value="true">검출 (양성)</option>
                <option value="false">미검출 (음성)</option>
                <option value="">판정 보류</option>
              </select>
            </label>
            <label>Ct 값 (선택)
              <input id="grpVfCt" type="text" placeholder="예: 23.4">
            </label>
            <label>형광물질(채널)
              <input id="grpVfFluorophore" type="text" placeholder="예: FAM">
            </label>
          </div>
          <label style="display:block;margin-top:10px;">메모
            <textarea id="grpVfNote" rows="2" placeholder="실험 조건, 특이사항 등" style="width:100%;resize:vertical;padding:8px 10px;border-radius:8px;border:1px solid var(--border);background:var(--bg);color:var(--fg);font-size:12.5px;"></textarea>
          </label>
          <div class="form-actions" style="margin-top:12px;">
            <button id="grpVfSaveBtn" class="btn btn-primary">검증 결과 저장</button>
            <span id="grpVfSaveStatus" class="form-status"></span>
          </div>
        </div>
      </div>
    </div>
  `;
  detail.querySelectorAll(".group-chip").forEach((chip) => {
    chip.addEventListener("click", () => selectSet(chip.dataset.setId));
  });
  const treeCard = detail.querySelector(".tree-card");
  $("#loadTreeBtn").addEventListener("click", () => renderTreeSection($("#treeContainer"), g, g.segment));
  $("#treeFullscreenBtn").addEventListener("click", () => toggleFullscreen(treeCard));
  $("#deleteGroupBtn").addEventListener("click", () => deleteGroup(g.group_id, g.label));
  $("#renameGroupBtn").addEventListener("click", () => renameGroup(g.group_id, g.label));

  $("#saveNoteBtn").addEventListener("click", async () => {
    const statusEl = $("#noteSaveStatus");
    statusEl.textContent = "저장 중...";
    try {
      const note = $("#setNoteInput").value;
      const res = await fetch(`/api/groups/${encodeURIComponent(g.group_id)}/note`, {
        method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ note }),
      });
      if (!res.ok) throw new Error(`HTTP ${res.status}`);
      statusEl.textContent = "저장됨";
    } catch (e) {
      statusEl.textContent = `저장 실패: ${e.message}`;
    }
  });
  wireNoteImages(g, "groups");
  wireGroupVerification(g);
  loadSetDiversityOverlay(g);

  const setsById = {};
  (g.member_details || []).forEach((s) => { setsById[s.set_id] = s; });
  wireComboSection("combo", g.combo_genomes || [], setsById);
  wireComboSection("standardCombo", g.standard_combo || [], setsById);
  wireComboSection("refgenomeCombo", g.refgenome_combo || [], setsById);
}

async function deleteGroup(groupId, label) {
  if (!await uiConfirm(`그룹 "${label}"을(를) 삭제할까요? (구성 세트 자체는 삭제되지 않습니다)`)) return;
  try {
    const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}`, { method: "DELETE" });
    if (!res.ok) {
      const data = await res.json().catch(() => ({}));
      throw new Error(data.error || `HTTP ${res.status}`);
    }
    await loadGroups();
    showOverview();
  } catch (e) {
    uiAlert(`그룹 삭제 실패: ${e.message}`);
  }
}

async function renameGroup(groupId, currentLabel) {
  const name = await uiPrompt("새 그룹 이름을 입력하세요:", currentLabel);
  if (name === null || !name.trim() || name.trim() === currentLabel) return;
  try {
    const res = await fetch(`/api/groups/${encodeURIComponent(groupId)}/rename`, {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name: name.trim() }),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    const g = state.groups.find((x) => x.group_id === groupId);
    if (g) g.label = data.label;
    renderSidebar();
    if (state.activeType === "group" && state.activeId === groupId) {
      renderGroupDetail(await fetchJSON(`/api/groups/${encodeURIComponent(groupId)}`));
    } else if (!$("#overview").hidden) {
      showOverview();
    }
  } catch (e) {
    uiAlert(`이름 변경 실패: ${e.message}`);
  }
}

async function renameSet(setId, currentName) {
  const name = await uiPrompt("새 이름을 입력하세요 (원래 ID로 되돌리려면 비워두고 확인):", currentName);
  if (name === null) return;
  try {
    const res = await fetch(`/api/sets/${encodeURIComponent(setId)}/rename`, {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ name: name.trim() }),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    const s = state.sets.find((x) => x.set_id === setId);
    if (s) s.display_name = data.display_name;
    renderSidebar();
    if (state.activeType === "set" && state.activeId === setId) {
      renderDetail(await fetchJSON(`/api/sets/${encodeURIComponent(setId)}`));
    } else if (!$("#overview").hidden) {
      showOverview();
    }
  } catch (e) {
    uiAlert(`이름 변경 실패: ${e.message}`);
  }
}

// ---------------- add-set form ----------------

function showAddForm() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("addForm");
  $("#fStatus").textContent = "";
}

async function showParamsPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("paramsView");
  const statusEl = $("#paramsSaveStatus");
  statusEl.textContent = "";
  try {
    const p = await fetchJSON("/api/params");
    $("#pTargetCoverage").value = p.target_coverage_pct;
    $("#pPrimerMismatch").value = p.primer_mismatch_pass;
    $("#pProbeMismatch").value = p.probe_mismatch;
    $("#pThreePrime").value = p.three_prime_consecutive_match;
    $("#pMfeprimerCpu").value = p.mfeprimer_cpu;
    $("#pMfeprimerBin").value = p.mfeprimer_bin || "";
    showMfeprimerStatus(p.mfeprimer_status, !p.mfeprimer_bin);
    const t = p.tools || {};
    $("#toolsInfo").innerHTML = `conda 환경: <code>${escHtml(t.env_prefix || "(미지정 — PATH 사용)")}</code> &middot; snakemake: <code>${escHtml(t.snakemake || "찾을 수 없음")}</code> &middot; mafft: <code>${escHtml(t.mafft || "찾을 수 없음")}</code>`;
  } catch (e) {
    statusEl.textContent = `불러오기 실패: ${e.message}`;
  }
}

function showMfeprimerStatus(st, auto) {
  const el = $("#mfeprimerStatus");
  if (!st) { el.textContent = ""; return; }
  el.style.color = st.ok ? "var(--ok)" : "var(--bad)";
  el.innerHTML = st.ok
    ? `✔ ${escHtml(st.version)} — <code>${escHtml(st.path)}</code>${auto ? " (자동 탐색)" : ""}`
    : `✘ ${escHtml(st.error || "mfeprimer를 찾을 수 없습니다")}${st.path ? ` — <code>${escHtml(st.path)}</code>` : ""}`;
}

$("#checkMfeprimerBtn").addEventListener("click", async () => {
  const path = $("#pMfeprimerBin").value.trim();
  $("#mfeprimerStatus").textContent = "확인 중...";
  try {
    showMfeprimerStatus(await fetchJSON(`/api/params/check-mfeprimer?path=${encodeURIComponent(path)}`));
  } catch (e) {
    $("#mfeprimerStatus").textContent = `확인 실패: ${e.message}`;
  }
});

$("#saveParamsBtn").addEventListener("click", async () => {
  const statusEl = $("#paramsSaveStatus");
  statusEl.textContent = "저장 중...";
  try {
    const payload = {
      target_coverage_pct: Number($("#pTargetCoverage").value),
      primer_mismatch_pass: Number($("#pPrimerMismatch").value),
      probe_mismatch: Number($("#pProbeMismatch").value),
      three_prime_consecutive_match: Number($("#pThreePrime").value),
      mfeprimer_cpu: Number($("#pMfeprimerCpu").value),
      mfeprimer_bin: $("#pMfeprimerBin").value.trim(),
    };
    const res = await fetch("/api/params", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    showMfeprimerStatus(data.mfeprimer_status, !data.mfeprimer_bin);
    statusEl.textContent = "저장됨. 기존 세트에 반영하려면 \"분석 다시하기\"를 눌러주세요.";
  } catch (e) {
    statusEl.textContent = `저장 실패: ${e.message}`;
  }
});

$("#paramsBtn").addEventListener("click", showParamsPage);

// ---------------- changelog page ----------------

function renderMarkdownLite(md) {
  // Just enough markdown for CHANGELOG.md's own syntax (headers, "## [x.y.z] - date"
  // version headers, bullet lists, `code`, **bold**, [text](url) links) -- not a general
  // markdown renderer, so it never has to handle anything CHANGELOG.md itself doesn't use.
  const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
  const inline = (s) => esc(s)
    .replace(/`([^`]+)`/g, "<code>$1</code>")
    .replace(/\*\*([^*]+)\*\*/g, "<strong>$1</strong>")
    .replace(/\[([^\]]+)\]\(([^)]+)\)/g, '<a href="$2" target="_blank" rel="noopener">$1</a>');
  const lines = md.split("\n");
  let html = "";
  let inList = false;
  const closeList = () => { if (inList) { html += "</ul>"; inList = false; } };
  for (const raw of lines) {
    const line = raw.trimEnd();
    if (/^### /.test(line)) { closeList(); html += `<h3>${inline(line.slice(4))}</h3>`; }
    else if (/^## /.test(line)) { closeList(); html += `<h2 class="changelog-version">${inline(line.slice(3))}</h2>`; }
    else if (/^# /.test(line)) { closeList(); html += `<h1>${inline(line.slice(2))}</h1>`; }
    else if (/^-\s+/.test(line)) {
      if (!inList) { html += "<ul>"; inList = true; }
      html += `<li>${inline(line.replace(/^-\s+/, ""))}</li>`;
    } else if (line.trim() === "") {
      closeList();
    } else {
      closeList();
      html += `<p>${inline(line)}</p>`;
    }
  }
  closeList();
  return html;
}

async function showChangelogPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("changelogView");
  const view = $("#changelogView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    const { markdown } = await fetchJSON("/api/changelog");
    view.innerHTML = `<div class="changelog-body">${renderMarkdownLite(markdown)}</div>`;
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
  }
}

$("#changelogBtn").addEventListener("click", showChangelogPage);

// ---------------------------------------------------------------------------------------
// XlGrid: an Excel-style grid ("스프레드시트 형식") shared by the primer and standard-material pages.
// It edits the SAME row objects the classic form-style table uses (opts.getRows()), so saving, validation
// and the server round trip are unchanged -- the button only switches how the rows are shown and edited.
//   cells      : click / drag / Shift+click / Shift+arrows select a range; row numbers and column letters
//                select whole rows / columns; Ctrl+A selects everything
//   editing    : type to replace, F2 or double-click to edit, Enter/Tab commit and move, Esc cancels;
//                the bar above the grid shows the active cell's full content and edits it too (long sequences)
//   clipboard  : Ctrl+C / X / V as tab-separated text (Excel, Google Sheets); one value pasted over a range fills it
//   other      : Delete clears, Ctrl+D fills down, Ctrl+Z / Ctrl+Y undo / redo, "True/False" cells toggle with Space
const XL_ROW_NO_W = 46;
const xlLetters = (n) => { let s = ""; n += 1; while (n > 0) { const m = (n - 1) % 26; s = String.fromCharCode(65 + m) + s; n = Math.floor((n - 1) / 26); } return s; };
const XL_TRUE = (v) => ["true", "1", "yes", "y"].includes(String(v).trim().toLowerCase());

function createXlGrid(opts) {
  const host = opts.host, cols = opts.columns;
  const rows = () => opts.getRows();
  let anchor = { r: 0, c: 0 }, focus = { r: 0, c: 0 }, dragging = false, editor = null, hist = [], redo = [];

  const colValue = (col, row) => (col.value ? col.value(row) : (row[col.key] ?? ""));
  const colText = (col, row) => { const v = colValue(col, row); return col.type === "bool" ? (XL_TRUE(v) ? "True" : "False") : String(v); };
  const normalize = (col, v) => {
    v = String(v ?? "");
    if (col.type === "bool") return XL_TRUE(v) || v.trim() === "" ? "True" : "False";
    if (col.seq) return v.replace(/\s+/g, "").toUpperCase();
    if (col.upper) return v.trim().toUpperCase();
    return v.replace(/[\t\r\n]+/g, " ").trim();
  };
  const range = () => ({ r0: Math.min(anchor.r, focus.r), r1: Math.max(anchor.r, focus.r), c0: Math.min(anchor.c, focus.c), c1: Math.max(anchor.c, focus.c) });
  const clampR = (r) => Math.max(0, Math.min(rows().length - 1, r));
  const clampC = (c) => Math.max(0, Math.min(cols.length - 1, c));
  const stickyLeft = []; { let x = XL_ROW_NO_W; cols.forEach((c, i) => { stickyLeft[i] = c.sticky ? x : null; if (c.sticky) x += c.width; }); }

  host.innerHTML = `<div class="xl-bar"><span class="xl-addr">A1</span><input class="xl-formula" type="text" spellcheck="false" autocomplete="off"><span class="xl-err"></span></div>
    <div class="xl-tools"><button type="button" class="btn btn-secondary" data-xl="undo" title="실행 취소 (Ctrl+Z)">&#8630; 실행 취소</button>
      <button type="button" class="btn btn-secondary" data-xl="redo" title="다시 실행 (Ctrl+Y)">&#8631; 다시 실행</button>
      <button type="button" class="btn btn-secondary" data-xl="insert" title="선택한 행 아래에 빈 행 삽입">+ 행 삽입</button>
      <button type="button" class="btn btn-secondary" data-xl="delete" title="선택한 행 삭제">&#8722; 선택 행 삭제</button></div>
    <div class="xl-wrap" tabindex="0"></div>`;
  const wrap = host.querySelector(".xl-wrap"), bar = host.querySelector(".xl-formula"), addr = host.querySelector(".xl-addr"), errEl = host.querySelector(".xl-err");

  function renderBody() {
    const errs = opts.getErrors ? opts.getErrors() : {};
    const head1 = `<tr class="xl-letters"><th class="xl-corner" rowspan="2" data-corner></th>${cols.map((c, i) => `<th class="xl-ch${stickyLeft[i] !== null ? " xl-st" : ""}" data-ch="${i}" ${stickyLeft[i] !== null ? `style="left:${stickyLeft[i]}px"` : ""}>${xlLetters(i)}</th>`).join("")}</tr>`;
    const head2 = `<tr class="xl-titles">${cols.map((c, i) => `<th class="xl-title${stickyLeft[i] !== null ? " xl-st" : ""}" data-ch="${i}" ${stickyLeft[i] !== null ? `style="left:${stickyLeft[i]}px"` : ""}>${escHtml(c.title)}</th>`).join("")}</tr>`;
    const body = rows().map((row, r) => `<tr data-r="${r}"><th class="xl-rh" data-rh="${r}">${r + 1}</th>${cols.map((c, i) => {
      const txt = colText(c, row), e = errs[r] && errs[r][c.key];
      const cls = ["xl-c", c.mono || c.seq ? "mono" : "", c.ro ? "ro" : "", c.type === "bool" ? (txt === "True" ? "xl-true" : "xl-false") : "", c.num ? "num" : "", stickyLeft[i] !== null ? "xl-st" : "", e ? "xl-bad" : ""].filter(Boolean).join(" ");
      return `<td class="${cls}" data-r="${r}" data-c="${i}" ${stickyLeft[i] !== null ? `style="left:${stickyLeft[i]}px"` : ""} ${e ? `title="${escHtml(e)}"` : ""}>${escHtml(txt)}</td>`;
    }).join("")}</tr>`).join("");
    wrap.innerHTML = `<table class="xl-grid"><colgroup><col style="width:${XL_ROW_NO_W}px">${cols.map((c) => `<col style="width:${c.width}px">`).join("")}</colgroup><thead>${head1}${head2}</thead><tbody>${body}</tbody></table>`;
    paint();
  }

  function cellEl(r, c) { return wrap.querySelector(`td[data-r="${r}"][data-c="${c}"]`); }

  function paint() {
    if (!rows().length) return;
    anchor = { r: clampR(anchor.r), c: clampC(anchor.c) }; focus = { r: clampR(focus.r), c: clampC(focus.c) };
    const { r0, r1, c0, c1 } = range();
    wrap.querySelectorAll("td.sel, td.active, th.hl").forEach((e) => e.classList.remove("sel", "active", "hl"));
    for (let r = r0; r <= r1; r++) for (let c = c0; c <= c1; c++) { const td = cellEl(r, c); if (td) td.classList.add("sel"); }
    const act = cellEl(focus.r, focus.c);
    if (act) act.classList.add("active");
    for (let r = r0; r <= r1; r++) wrap.querySelector(`th[data-rh="${r}"]`)?.classList.add("hl");
    for (let c = c0; c <= c1; c++) wrap.querySelectorAll(`th[data-ch="${c}"]`).forEach((e) => e.classList.add("hl"));
    const col = cols[focus.c], row = rows()[focus.r];
    addr.textContent = `${xlLetters(focus.c)}${focus.r + 1}`;
    bar.value = colText(col, row);
    bar.readOnly = !!col.ro;
    bar.classList.toggle("ro", !!col.ro);
    const e = opts.getErrors && (opts.getErrors()[focus.r] || {})[col.key];
    errEl.textContent = e || "";
    host.querySelector('[data-xl="undo"]').disabled = !hist.length;
    host.querySelector('[data-xl="redo"]').disabled = !redo.length;
    if (act && !editor) {
      const w = wrap.getBoundingClientRect(), a = act.getBoundingClientRect();
      const topPad = 48, leftPad = XL_ROW_NO_W + (stickyLeft[focus.c] !== null ? 0 : cols.reduce((s, c2, i) => s + (stickyLeft[i] !== null ? c2.width : 0), 0));
      if (a.top < w.top + topPad) wrap.scrollTop -= (w.top + topPad - a.top);
      else if (a.bottom > w.bottom - 16) wrap.scrollTop += a.bottom - (w.bottom - 16);
      if (stickyLeft[focus.c] === null) {
        if (a.left < w.left + leftPad) wrap.scrollLeft -= (w.left + leftPad - a.left);
        else if (a.right > w.right - 16) wrap.scrollLeft += a.right - (w.right - 16);
      }
    }
  }

  const snapshot = () => { hist.push(JSON.stringify(rows())); if (hist.length > 100) hist.shift(); redo = []; };
  function restore(json) { const arr = rows(); arr.splice(0, arr.length, ...JSON.parse(json)); changed(); }
  function undo() { if (!hist.length) return; redo.push(JSON.stringify(rows())); restore(hist.pop()); }
  function redoIt() { if (!redo.length) return; hist.push(JSON.stringify(rows())); restore(redo.pop()); }
  function changed() { if (opts.ensureBlank) opts.ensureBlank(); if (opts.onChange) opts.onChange(); renderBody(); }

  function setCells(list, snapped) { // [{r, c, v}] -> one undo step (snapped: the caller already took the snapshot)
    let any = !!snapped;
    for (const { r, c, v } of list) {
      const col = cols[c]; if (col.ro || r >= rows().length) continue;
      const nv = normalize(col, v);
      if (String(rows()[r][col.key] ?? "") !== nv) { if (!any) { snapshot(); any = true; } rows()[r][col.key] = nv; }
    }
    if (any) changed(); else paint();
  }

  function startEdit(initial) {
    const col = cols[focus.c];
    if (col.ro) return;
    if (col.type === "bool") { setCells([{ r: focus.r, c: focus.c, v: XL_TRUE(rows()[focus.r][col.key]) ? "False" : "True" }]); return; }
    const td = cellEl(focus.r, focus.c);
    if (!td) return;
    const cur = String(rows()[focus.r][col.key] ?? "");
    td.classList.add("editing");
    td.innerHTML = `<input class="xl-editor${col.seq || col.mono ? " mono" : ""}" type="text" spellcheck="false" autocomplete="off">`;
    const inp = td.firstChild;
    inp.value = initial === undefined ? cur : initial;
    editor = { r: focus.r, c: focus.c, orig: cur, el: inp, done: false };
    inp.focus();
    if (initial === undefined) inp.select(); else inp.setSelectionRange(inp.value.length, inp.value.length);
    inp.addEventListener("input", () => { bar.value = inp.value; });
    inp.addEventListener("keydown", (ev) => {
      if (ev.key === "Enter") { ev.preventDefault(); finishEdit(true, ev.shiftKey ? -1 : 1, 0); }
      else if (ev.key === "Tab") { ev.preventDefault(); finishEdit(true, 0, ev.shiftKey ? -1 : 1); }
      else if (ev.key === "Escape") { ev.preventDefault(); finishEdit(false, 0, 0); }
      ev.stopPropagation();
    });
    inp.addEventListener("blur", () => { if (editor && editor.el === inp) finishEdit(true, 0, 0, true); });
  }

  function finishEdit(commit, dr, dc, fromBlur) {
    if (!editor || editor.done) return;
    const ed = editor; ed.done = true; editor = null;
    const v = ed.el.value;
    if (commit) setCells([{ r: ed.r, c: ed.c, v }]); else renderBody();
    if (dr || dc) move(dr, dc, false);
    if (!fromBlur) wrap.focus();
  }

  function move(dr, dc, extend) {
    focus = { r: clampR(focus.r + dr), c: clampC(focus.c + dc) };
    if (!extend) anchor = { ...focus };
    paint();
  }

  const toTsv = () => { const { r0, r1, c0, c1 } = range(); const out = []; for (let r = r0; r <= r1; r++) { const line = []; for (let c = c0; c <= c1; c++) line.push(colText(cols[c], rows()[r]).replace(/[\t\r\n]+/g, " ")); out.push(line.join("\t")); } return out.join("\n"); };
  function clearRange() { const { r0, r1, c0, c1 } = range(); const l = []; for (let r = r0; r <= r1; r++) for (let c = c0; c <= c1; c++) l.push({ r, c, v: cols[c].type === "bool" ? "True" : "" }); setCells(l); }
  function pasteText(text) {
    const lines = text.replace(/\r/g, "").replace(/\n$/, "").split("\n").map((l) => l.split("\t"));
    if (!lines.length || (lines.length === 1 && lines[0].length === 1 && lines[0][0] === "" )) return;
    const { r0, r1, c0, c1 } = range(), list = [];
    let snapped = false;
    if (lines.length === 1 && lines[0].length === 1 && (r1 > r0 || c1 > c0)) {
      for (let r = r0; r <= r1; r++) for (let c = c0; c <= c1; c++) list.push({ r, c, v: lines[0][0] });
    } else {
      lines.forEach((cells, li) => cells.forEach((v, ci) => { if (c0 + ci < cols.length) list.push({ r: r0 + li, c: c0 + ci, v }); }));
      const need = r0 + lines.length;
      if (need > rows().length) { snapshot(); snapped = true; while (rows().length < need) rows().push(opts.blankRow()); }
      anchor = { r: r0, c: c0 }; focus = { r: Math.min(need - 1, rows().length - 1), c: Math.min(cols.length - 1, c0 + Math.max(...lines.map((l) => l.length)) - 1) };
    }
    setCells(list, snapped);
  }

  // ---- events
  wrap.addEventListener("mousedown", (ev) => {
    if (ev.target.closest(".xl-editor")) return;
    if (editor) finishEdit(true, 0, 0, true);
    const td = ev.target.closest("td.xl-c"), rh = ev.target.closest("th.xl-rh"), ch = ev.target.closest("th.xl-ch, th.xl-title"), corner = ev.target.closest("th.xl-corner");
    if (td) { focus = { r: +td.dataset.r, c: +td.dataset.c }; if (!ev.shiftKey) anchor = { ...focus }; dragging = "cell"; }
    else if (rh) { const r = +rh.dataset.rh; if (!ev.shiftKey) anchor = { r, c: 0 }; focus = { r, c: cols.length - 1 }; dragging = "row"; }
    else if (ch) { const c = +ch.dataset.ch; if (!ev.shiftKey) anchor = { r: 0, c }; focus = { r: rows().length - 1, c }; dragging = "col"; }
    else if (corner) { anchor = { r: 0, c: 0 }; focus = { r: rows().length - 1, c: cols.length - 1 }; }
    else return;
    paint(); wrap.focus(); ev.preventDefault();
  });
  wrap.addEventListener("mouseover", (ev) => {
    if (!dragging) return;
    const td = ev.target.closest("td.xl-c"), rh = ev.target.closest("th.xl-rh"), ch = ev.target.closest("th.xl-ch, th.xl-title");
    if (dragging === "cell" && td) focus = { r: +td.dataset.r, c: +td.dataset.c };
    else if (dragging === "row" && (rh || td)) focus = { r: +(rh ? rh.dataset.rh : td.dataset.r), c: cols.length - 1 };
    else if (dragging === "col" && (ch || td)) focus = { r: rows().length - 1, c: +(ch ? ch.dataset.ch : td.dataset.c) };
    else return;
    paint();
  });
  const stopDrag = () => { dragging = false; };
  document.addEventListener("mouseup", stopDrag);
  wrap.addEventListener("dblclick", (ev) => { if (ev.target.closest("td.xl-c")) startEdit(); });

  wrap.addEventListener("keydown", (ev) => {
    if (editor) return;
    const mod = ev.ctrlKey || ev.metaKey, k = ev.key;
    if (mod && !ev.altKey) {
      const lk = k.toLowerCase();
      if (lk === "a") { ev.preventDefault(); anchor = { r: 0, c: 0 }; focus = { r: rows().length - 1, c: cols.length - 1 }; paint(); return; }
      if (lk === "z") { ev.preventDefault(); ev.shiftKey ? redoIt() : undo(); return; }
      if (lk === "y") { ev.preventDefault(); redoIt(); return; }
      if (lk === "d") { ev.preventDefault(); const { r0, r1, c0, c1 } = range(), l = []; for (let r = r0 + 1; r <= r1; r++) for (let c = c0; c <= c1; c++) l.push({ r, c, v: colText(cols[c], rows()[r0]) }); setCells(l); return; }
      if (["c", "x", "v"].includes(lk)) return;   // handled by the clipboard events below
      if (k.startsWith("Arrow")) { ev.preventDefault(); const to = { ArrowUp: [-rows().length, 0], ArrowDown: [rows().length, 0], ArrowLeft: [0, -cols.length], ArrowRight: [0, cols.length] }[k]; move(to[0], to[1], ev.shiftKey); return; }
      return;
    }
    const arrows = { ArrowUp: [-1, 0], ArrowDown: [1, 0], ArrowLeft: [0, -1], ArrowRight: [0, 1] };
    if (arrows[k]) { ev.preventDefault(); move(arrows[k][0], arrows[k][1], ev.shiftKey); }
    else if (k === "Tab") { ev.preventDefault(); const c = focus.c + (ev.shiftKey ? -1 : 1); if (c >= cols.length) { focus = { r: clampR(focus.r + 1), c: 0 }; anchor = { ...focus }; paint(); } else if (c < 0) { focus = { r: clampR(focus.r - 1), c: cols.length - 1 }; anchor = { ...focus }; paint(); } else move(0, ev.shiftKey ? -1 : 1, false); }
    else if (k === "Enter") { ev.preventDefault(); move(ev.shiftKey ? -1 : 1, 0, false); }
    else if (k === "F2") { ev.preventDefault(); startEdit(); }
    else if (k === "Home") { ev.preventDefault(); move(0, -cols.length, ev.shiftKey); }
    else if (k === "End") { ev.preventDefault(); move(0, cols.length, ev.shiftKey); }
    else if (k === "PageDown") { ev.preventDefault(); move(10, 0, ev.shiftKey); }
    else if (k === "PageUp") { ev.preventDefault(); move(-10, 0, ev.shiftKey); }
    else if (k === "Delete" || k === "Backspace") { ev.preventDefault(); clearRange(); }
    else if (k === " " && cols[focus.c].type === "bool") { ev.preventDefault(); startEdit(); }
    else if (k.length === 1 && !ev.altKey) { ev.preventDefault(); startEdit(k); }
  });
  wrap.addEventListener("copy", (ev) => { if (editor) return; ev.clipboardData.setData("text/plain", toTsv()); ev.preventDefault(); });
  wrap.addEventListener("cut", (ev) => { if (editor) return; ev.clipboardData.setData("text/plain", toTsv()); ev.preventDefault(); clearRange(); });
  wrap.addEventListener("paste", (ev) => { if (editor) return; ev.preventDefault(); pasteText((ev.clipboardData || window.clipboardData).getData("text")); });

  bar.addEventListener("keydown", (ev) => {
    if (ev.key === "Enter") { ev.preventDefault(); if (!bar.readOnly) setCells([{ r: focus.r, c: focus.c, v: bar.value }]); wrap.focus(); }
    else if (ev.key === "Escape") { ev.preventDefault(); paint(); wrap.focus(); }
  });
  host.querySelector(".xl-tools").addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-xl]"); if (!b) return;
    const act = b.dataset.xl;
    if (act === "undo") undo(); else if (act === "redo") redoIt();
    else if (act === "insert") { snapshot(); const at = range().r1 + 1; rows().splice(at, 0, opts.blankRow()); anchor = focus = { r: at, c: focus.c }; anchor = { ...focus }; changed(); }
    else if (act === "delete") {
      const { r0, r1 } = range(); snapshot();
      rows().splice(r0, r1 - r0 + 1);
      if (!rows().length) rows().push(opts.blankRow());
      focus = { r: clampR(r0), c: focus.c }; anchor = { ...focus }; changed();
    }
    wrap.focus();
  });

  renderBody();
  return { refresh: renderBody, focus: () => wrap.focus(), destroy: () => document.removeEventListener("mouseup", stopDrag) };
}

// Per-page choice of grid vs classic view, remembered in this browser.
const XL_MODE_KEY = "lasv_xl_mode_v1";
function xlModeOn(page) { try { return !!JSON.parse(localStorage.getItem(XL_MODE_KEY) || "{}")[page]; } catch { return false; } }
function xlModeSet(page, on) { try { const m = JSON.parse(localStorage.getItem(XL_MODE_KEY) || "{}"); m[page] = !!on; localStorage.setItem(XL_MODE_KEY, JSON.stringify(m)); } catch { /* private mode */ } }

// ---------------- primers manifest spreadsheet (config/primers.tsv) ----------------
// A lightweight, dependency-free editable grid -- no spreadsheet library anywhere else in
// this codebase, and the row count (tens, not thousands) doesn't need one. The whole table
// is held in `primersState.rows` and re-rendered on every structural change (add/delete
// row); per-cell edits just write back into that in-memory row object on blur/change, no
// re-render needed until save.

let primersState = { rows: [], errorsByRow: {} };

function primersBlankRow() {
  return { primer_id: "", primer_name: "", segment: "S", forward_primer: "", reverse_primer: "", probe: "",
           fluorophore: "", source: "", remark: "", proceed: "True", add_date: "" };
}

async function showPrimersPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("primersView");
  const view = $("#primersView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    const data = await fetchJSON("/api/primers/manifest");
    primersState = { rows: data.rows || [], errorsByRow: {}, savedIds: new Map((data.rows || []).filter((r) => r.primer_id).map((r) => [r.primer_id, r.primer_name])) };
    // Always keep one trailing blank row ready to type a new entry into.
    if (!primersState.rows.length || primersState.rows[primersState.rows.length - 1].primer_name) {
      primersState.rows.push(primersBlankRow());
    }
    renderPrimersPage();
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
  }
}

const PRIMERS_FIELDS = ["primer_name", "segment", "forward_primer", "reverse_primer", "probe", "fluorophore", "source", "remark", "proceed"];

function primersRowHtml(r, i) {
  const errs = primersState.errorsByRow[i] || {};
  const ec = (f) => (errs[f] ? " cell-error" : "");
  const proceedOn = ["true", "1", "yes", "y"].includes(String(r.proceed).toLowerCase());
  const txt = (f, extra = "", ph = "") =>
    `<td class="${ec(f).trim()}"><input type="text" class="cell-input ${extra}" data-field="${f}" value="${escHtml(r[f] || "")}" placeholder="${ph}" spellcheck="false" autocomplete="off"></td>`;
  const errRow = Object.keys(errs).length
    ? `<tr class="err-row"><td class="row-num"></td><td colspan="12"><div class="primers-cell-errors">${Object.values(errs).map(escHtml).join(" / ")}</div></td></tr>`
    : "";
  return `
    <tr data-row="${i}" class="${r.add_date ? "" : "row-new"}">
      <td class="row-num">${i + 1}</td>
      <td class="primer-id-cell mono" title="영구 ID — 이름을 바꿔도 변하지 않고, 결과 폴더(output/primers/ID)와 메모·검증 기록의 키입니다">${escHtml(r.primer_id || "자동")}</td>
      ${txt("primer_name", "", "예: MyLab_S_Test1")}
      <td class="${ec("segment").trim()}"><select class="cell-input" data-field="segment">
        <option value="S" ${r.segment === "S" ? "selected" : ""}>S</option>
        <option value="L" ${r.segment === "L" ? "selected" : ""}>L</option>
      </select></td>
      ${txt("forward_primer", "seq-cell", "")}
      ${txt("reverse_primer", "seq-cell", "")}
      ${txt("probe", "seq-cell", "")}
      ${txt("fluorophore", "", "FAM")}
      ${txt("source", "", "")}
      ${txt("remark", "remark-cell", "제조사 제품명 등")}
      <td class="center-cell"><input type="checkbox" class="cell-input" data-field="proceed" ${proceedOn ? "checked" : ""} title="해제하면 이 세트는 분석에서 완전히 제외됩니다"></td>
      <td class="add-date-cell">${escHtml(r.add_date || "저장 시 자동")}</td>
      <td class="del-cell"><button class="row-del-btn" data-del-row="${i}" title="행 삭제">&times;</button></td>
    </tr>${errRow}`;
}

function primersFocusCell(rowIdx, field) {
  const el = document.querySelector(`#primersView tr[data-row="${rowIdx}"] [data-field="${field}"]`);
  if (el) { el.focus(); if (el.select) el.select(); }
}

function primersEnsureTrailingBlank() {
  const last = primersState.rows[primersState.rows.length - 1];
  if (!last || last.primer_name || last.forward_primer || last.reverse_primer || last.probe) {
    primersState.rows.push(primersBlankRow());
    return true;
  }
  return false;
}

function primersXlColumns() {
  return [
    { key: "primer_id", title: "ID", width: 66, ro: true, mono: true, sticky: true, value: (r) => r.primer_id || "자동" },
    { key: "primer_name", title: "이름 (primer_name)", width: 240, sticky: true },
    { key: "segment", title: "Seg", width: 52, upper: true },
    { key: "forward_primer", title: "Forward Primer", width: 230, seq: true },
    { key: "reverse_primer", title: "Reverse Primer", width: 230, seq: true },
    { key: "probe", title: "Probe", width: 230, seq: true },
    { key: "fluorophore", title: "형광물질", width: 100 },
    { key: "source", title: "출처", width: 110 },
    { key: "remark", title: "비고", width: 190 },
    { key: "proceed", title: "분석진행", width: 84, type: "bool" },
    { key: "add_date", title: "추가일", width: 100, ro: true, value: (r) => r.add_date || "저장 시 자동" },
  ];
}

// Excel (.xlsx) import: the file is parsed on the server and opened in a spreadsheet pop-up where every cell can be corrected;
// the same checks as saving run live (errors are marked in the cells) and "표에 적용" is possible once none is left.
// Nothing is saved until the page's own "저장" is pressed. `ctx` = { existing(), columns(), blankRow(), apply(rows, mode) }.
async function importXlsxInto(kind, file, ctx) {
  let data;
  try {
    const res = await fetch(`/api/xlsx/parse?kind=${kind}`, { method: "POST", headers: { "Content-Type": "application/octet-stream" }, body: await file.arrayBuffer() });
    data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
  } catch (e) { await uiAlert(`엑셀 파일을 읽지 못했습니다: ${e.message}`); return; }
  if (!data.rows.length) { await uiAlert("엑셀에서 읽을 행이 없습니다." + (data.warnings.length ? `\n${data.warnings.join("\n")}` : "")); return; }
  return openXlsxReviewModal(kind, data, ctx);
}

function openXlsxReviewModal(kind, data, ctx) {
  const isBlank = (r) => !Object.entries(r).some(([k, v]) => v && !["proceed", "segment", "add_date", "primer_id"].includes(k));
  const st = { rows: [...data.rows, ctx.blankRow()], mode: "append", errors: {}, extra: [], timer: null, seq: 0 };
  const overlay = el("div", "modal-overlay detail-modal-overlay");
  overlay.innerHTML = `<div class="modal-box detail-modal" role="dialog" aria-modal="true" style="width:min(1250px,96vw);">
    <div class="modal-head"><div><div class="modal-title">엑셀 업로드 확인 · 수정</div>
      <div class="modal-sub">${data.rows.length}행을 읽었습니다. 셀을 눌러 바로 고치고, 오류(빨간 표시)가 없어지면 표에 적용할 수 있습니다.</div></div></div>
    <div class="detail-modal-body">
      ${data.warnings.length ? `<ul class="ui-dialog-list">${data.warnings.map((w) => `<li>${escHtml(w)}</li>`).join("")}</ul>` : ""}
      <div style="display:flex;gap:18px;align-items:center;flex-wrap:wrap;margin:6px 0;">
        <label><input type="radio" name="xlsxMode" value="append" checked> 현재 목록 아래에 <b>추가</b> (새 항목으로 등록)</label>
        <label><input type="radio" name="xlsxMode" value="replace"> 표 전체를 이 내용으로 <b>교체</b> (엑셀에 없는 기존 행은 저장 시 삭제)</label>
        <span id="xlsxStatus" class="form-status"></span>
      </div>
      <div id="xlsxExtra" class="primers-cell-errors" style="margin:6px 0;" hidden></div>
      <div id="xlsxGrid" class="xl-host"></div>
    </div>
    <div class="ui-dialog-actions"><button type="button" class="btn btn-secondary" id="xlsxCancel">취소</button>
      <button type="button" class="btn btn-primary" id="xlsxApply" disabled>표에 적용</button></div></div>`;
  modalHost().appendChild(overlay);
  const close = () => { clearTimeout(st.timer); overlay.remove(); };
  const payloadRows = () => {
    const mine = st.rows.filter((r) => !isBlank(r));
    const base = st.mode === "append" ? ctx.existing() : [];
    return { base, mine };
  };
  const refresh = async () => {
    const { base, mine } = payloadRows();
    const my = ++st.seq;
    let errors = [];
    try {
      const res = await fetch(`/api/xlsx/validate?kind=${kind}`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ rows: [...base, ...mine] }) });
      errors = (await res.json()).errors || [];
    } catch (e) { errors = [{ row: -1, field: "", message: `검증 실패: ${e.message}` }]; }
    if (my !== st.seq) return;                                   // a newer edit already triggered another check
    st.errors = {}; st.extra = [];
    const idxMap = new Map();                                    // index among non-blank rows -> index in the grid
    let k = 0; st.rows.forEach((r, gi) => { if (!isBlank(r)) idxMap.set(k++, gi); });
    for (const e of errors) {
      const gi = idxMap.get(e.row - base.length);
      if (e.row >= base.length && gi !== undefined) { st.errors[gi] = st.errors[gi] || {}; st.errors[gi][e.field] = e.message; }
      else st.extra.push(e.row < 0 ? e.message : `기존 행 ${e.row + 1}: ${e.message}`);
    }
    const n = errors.length;
    overlay.querySelector("#xlsxStatus").textContent = n ? `⚠ 오류 ${n}건 — 빨간 칸을 고쳐 주세요` : `✔ 문제 없음 (${mine.length}행)`;
    overlay.querySelector("#xlsxStatus").style.color = n ? "var(--bad)" : "var(--ok)";
    overlay.querySelector("#xlsxApply").disabled = n > 0 || mine.length === 0;
    const ex = overlay.querySelector("#xlsxExtra"); ex.hidden = !st.extra.length; ex.innerHTML = st.extra.map(escHtml).join("<br>");
    if (!overlay.querySelector(".xl-cell-editor, .xl-editor, .xl-wrap input")) grid.refresh();      // do not redraw under a cell being typed in
  };
  const schedule = () => { clearTimeout(st.timer); st.timer = setTimeout(refresh, 350); };
  const grid = createXlGrid({
    host: overlay.querySelector("#xlsxGrid"), columns: ctx.columns(), getRows: () => st.rows, blankRow: ctx.blankRow,
    getErrors: () => st.errors,
    ensureBlank: () => { const l = st.rows[st.rows.length - 1]; if (!l || !isBlank(l)) st.rows.push(ctx.blankRow()); },
    onChange: schedule,
  });
  overlay.querySelectorAll('input[name="xlsxMode"]').forEach((r) => r.addEventListener("change", () => { st.mode = r.value; refresh(); }));
  overlay.querySelector("#xlsxCancel").addEventListener("click", close);
  overlay.querySelector("#xlsxApply").addEventListener("click", () => { const mine = st.rows.filter((r) => !isBlank(r)); close(); ctx.apply(mine, st.mode); });
  refresh();
}
function wireXlsxButtons(kind, upBtn, fileInput, ctx) {
  upBtn.addEventListener("click", () => fileInput.click());
  fileInput.addEventListener("change", async () => {
    const f = fileInput.files[0];
    fileInput.value = "";
    if (f) await importXlsxInto(kind, f, ctx);
  });
}

function renderPrimersPage() {
  const view = $("#primersView");
  const xlOn = xlModeOn("primers");
  view.innerHTML = `
    <h2>프라이머 입력 / 관리</h2>
    <div class="note">
      엑셀처럼 셀을 클릭해 바로 입력하세요 — <b>Enter/↓/↑</b>로 위아래 이동, <b>Tab</b>으로 옆 칸 이동,
      엑셀/구글시트에서 복사한 영역을 <b>붙여넣기(Ctrl/Cmd+V)</b>하면 여러 행·열이 한 번에 채워집니다.
      <b>분석진행</b> 체크를 해제한 행은 파이프라인 분석에서 완전히 제외됩니다. <b>저장하는 순간 그 세트의 기존 결과(결과 폴더, 요약 표의 행, Probe 설계 기록)가 삭제되고 대시보드에서 숨겨집니다</b>(서열은 이 표에 남으므로 다시 체크하고 분석하면 복구됩니다).
      저장 후 실제 분석에 반영하려면 사이드바의 "분석 다시하기"를 눌러주세요.
    </div>
    <div class="primers-toolbar">
      <button id="primersAddRowBtn" class="btn btn-secondary">+ 행 추가</button>
      <button id="primersXlBtn" class="btn btn-secondary xl-mode-btn${xlOn ? " is-on" : ""}" title="엑셀처럼 셀을 선택·복사·붙여넣기·실행 취소할 수 있는 격자 화면으로 전환">스프레드시트 형식</button>
      <a id="primersXlsxDl" class="btn btn-secondary" href="/api/primers/xlsx" download title="현재 프라이머 목록을 엑셀(.xlsx)로 내려받습니다 (열 설명 시트 포함)">&#11015; 엑셀 다운로드</a>
      <button id="primersXlsxUp" class="btn btn-secondary" title="엑셀(.xlsx) 파일의 행을 이 표로 불러옵니다 (저장 전에 확인할 수 있습니다)">&#11014; 엑셀 업로드</button>
      <input type="file" id="primersXlsxFile" accept=".xlsx" hidden>
      <button id="primersSaveBtn" class="btn btn-primary">저장</button>
      <span id="primersSaveStatus" class="form-status"></span>
    </div>
    ${xlOn ? '<div id="primersXl" class="xl-host"></div>' : ""}
    <div class="primers-grid-wrap" ${xlOn ? "hidden" : ""}>
      <table class="primers-table">
        <thead><tr>
          <th class="row-num">#</th><th>ID</th><th>이름 (primer_name)</th><th>Seg</th><th>Forward Primer</th><th>Reverse Primer</th><th>Probe</th>
          <th>형광물질</th><th>출처</th><th>비고</th><th style="text-align:center;">분석진행</th><th>추가일</th><th></th>
        </tr></thead>
        <tbody id="primersTbody">${primersState.rows.map(primersRowHtml).join("")}</tbody>
      </table>
    </div>
    <div class="primers-hint">서열은 IUPAC 코드(A C G T R Y W S K M B D H V N)만 허용됩니다. 이름은 중복될 수 없고 segment는 S 또는 L입니다. ${xlOn ? "스프레드시트 형식: 셀 클릭·드래그로 선택, 타이핑/F2/더블클릭으로 편집, Ctrl+C/V 복사·붙여넣기, Delete 지우기, Ctrl+Z 실행 취소, 위쪽 수식 입력줄에서 긴 서열을 편집할 수 있습니다. " : ""}ID(P001…)는 저장 시 자동 부여되며 이름을 바꿔도 유지됩니다. 비고에는 제조사가 붙인 제품명/오더번호 등 자유 메모를 적으세요(분석에는 영향 없음).</div>`;

  $("#primersXlBtn").addEventListener("click", () => { xlModeSet("primers", !xlOn); renderPrimersPage(); });
  wireXlsxButtons("primers", $("#primersXlsxUp"), $("#primersXlsxFile"), {
    existing: () => primersState.rows.filter((r) => r.primer_name || r.forward_primer || r.reverse_primer || r.probe),
    columns: primersXlColumns, blankRow: primersBlankRow,
    apply: (rows, mode) => {
      const base = mode === "replace" ? [] : primersState.rows.filter((r) => r.primer_name || r.forward_primer || r.reverse_primer || r.probe);
      const add = rows.map((r) => (mode === "append" ? { ...r, primer_id: "", add_date: "" } : r));   // appended rows are new primers
      primersState.rows = [...base, ...add, primersBlankRow()];
      primersState.errorsByRow = {};
      renderPrimersPage();
      $("#primersSaveStatus").textContent = `엑셀에서 ${rows.length}행을 불러왔습니다 — 확인 후 저장을 눌러 반영하세요.`;
    },
  });
  if (xlOn) {
    $("#primersAddRowBtn").addEventListener("click", () => { primersState.rows.push(primersBlankRow()); renderPrimersPage(); });
    $("#primersSaveBtn").addEventListener("click", savePrimersManifest);
    createXlGrid({
      host: $("#primersXl"), columns: primersXlColumns(), getRows: () => primersState.rows, blankRow: primersBlankRow,
      getErrors: () => primersState.errorsByRow,
      ensureBlank: () => { const l = primersState.rows[primersState.rows.length - 1]; if (!l || l.primer_name || l.forward_primer || l.reverse_primer || l.probe) primersState.rows.push(primersBlankRow()); },
      onChange: () => { primersState.errorsByRow = {}; },
    });
    return;
  }

  const tbody = $("#primersTbody");

  // One set of delegated handlers on the tbody -- rows can be appended/removed without
  // rebinding anything or losing focus in the cell being typed into.
  tbody.addEventListener("change", (ev) => {
    const el = ev.target.closest("[data-field]");
    if (!el) return;
    const i = Number(el.closest("tr[data-row]").dataset.row);
    primersState.rows[i][el.dataset.field] = el.type === "checkbox" ? (el.checked ? "True" : "False") : el.value;
    if (i === primersState.rows.length - 1 && primersState.rows[i].primer_name) primersAppendBlankRow();
  });

  tbody.addEventListener("keydown", (ev) => {
    const el = ev.target.closest("[data-field]");
    if (!el) return;
    const i = Number(el.closest("tr[data-row]").dataset.row);
    const field = el.dataset.field;
    const isText = el.type === "text";
    let target = null;
    if (ev.key === "Enter" || (ev.key === "ArrowDown" && !(el.tagName === "SELECT"))) {
      target = ev.shiftKey && ev.key === "Enter" ? i - 1 : i + 1;
    } else if (ev.key === "ArrowUp" && el.tagName !== "SELECT") {
      target = i - 1;
    } else if (isText && ev.key === "ArrowRight" && el.selectionStart === el.value.length && el.selectionEnd === el.value.length) {
      const f = PRIMERS_FIELDS[PRIMERS_FIELDS.indexOf(field) + 1];
      if (f) { ev.preventDefault(); primersFocusCell(i, f); }
      return;
    } else if (isText && ev.key === "ArrowLeft" && el.selectionStart === 0 && el.selectionEnd === 0) {
      const f = PRIMERS_FIELDS[PRIMERS_FIELDS.indexOf(field) - 1];
      if (f) { ev.preventDefault(); primersFocusCell(i, f); }
      return;
    }
    if (target !== null) {
      ev.preventDefault();
      if (el.type !== "checkbox") el.dispatchEvent(new Event("change", { bubbles: true }));
      if (target >= primersState.rows.length) primersAppendBlankRow();
      if (target >= 0) primersFocusCell(target, field);
    }
  });

  // Paste a block copied from Excel/Google Sheets (tab-separated columns, newline-separated
  // rows) starting at the focused cell, in the same column order this grid uses.
  tbody.addEventListener("paste", (ev) => {
    const el = ev.target.closest("[data-field]");
    if (!el) return;
    const text = (ev.clipboardData || window.clipboardData).getData("text");
    if (!/[\t\n]/.test(text.replace(/\n$/, ""))) return;  // plain single value: let the browser paste normally
    ev.preventDefault();
    const startRow = Number(el.closest("tr[data-row]").dataset.row);
    const startCol = PRIMERS_FIELDS.indexOf(el.dataset.field);
    const lines = text.replace(/\r/g, "").replace(/\n$/, "").split("\n");
    lines.forEach((line, li) => {
      const ri = startRow + li;
      while (ri >= primersState.rows.length) primersState.rows.push(primersBlankRow());
      line.split("\t").forEach((cell, ci) => {
        const f = PRIMERS_FIELDS[startCol + ci];
        if (!f) return;
        let v = cell.trim();
        if (f === "segment") v = v.toUpperCase();
        if (f === "proceed") v = ["true", "1", "yes", "y"].includes(v.toLowerCase()) ? "True" : (v === "" ? "True" : "False");
        primersState.rows[ri][f] = v;
      });
    });
    primersEnsureTrailingBlank();
    primersState.errorsByRow = {};
    renderPrimersPage();
    $("#primersSaveStatus").textContent = `${lines.length}행을 붙여넣었습니다. 확인 후 저장을 눌러주세요.`;
  });

  tbody.addEventListener("click", (ev) => {
    const btn = ev.target.closest("[data-del-row]");
    if (!btn) return;
    primersState.rows.splice(Number(btn.dataset.delRow), 1);
    if (!primersState.rows.length) primersState.rows.push(primersBlankRow());
    primersState.errorsByRow = {};
    renderPrimersPage();
  });

  $("#primersAddRowBtn").addEventListener("click", () => {
    primersAppendBlankRow();
    primersFocusCell(primersState.rows.length - 1, "primer_name");
  });
  $("#primersSaveBtn").addEventListener("click", savePrimersManifest);
}

function primersAppendBlankRow() {
  const last = primersState.rows[primersState.rows.length - 1];
  if (last && !last.primer_name && !last.forward_primer && !last.reverse_primer && !last.probe) return;
  primersState.rows.push(primersBlankRow());
  const i = primersState.rows.length - 1;
  $("#primersTbody").insertAdjacentHTML("beforeend", primersRowHtml(primersState.rows[i], i));
}

async function savePrimersManifest() {
  // Renaming is safe now (results, notes, wet-lab records and group membership are keyed by the
  // permanent primer id) -- only deleting a row that already has an id orphans anything.
  const nowIds = new Set(primersState.rows.map((r) => r.primer_id).filter(Boolean));
  const lost = [...(primersState.savedIds || new Map())].filter(([id]) => !nowIds.has(id)).map(([id, n]) => `${id}  ${n}`);
  if (lost.length && !await uiConfirm(`다음 프라이머가 삭제됩니다:\n\n${lost.join("\n")}\n\n기존 분석 결과(output/primers/ID)·메모·실험 검증 기록·그룹 구성과의 연결이 끊어집니다. 계속 저장할까요?`)) {
    $("#primersSaveStatus").textContent = "저장 취소됨";
    return;
  }
  $("#primersSaveStatus").textContent = "저장 중...";
  toast("프라이머 표 저장 중...", { type: "run", id: "primers-save" });
  let res, data;
  try {
    res = await fetch("/api/primers/manifest", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ rows: primersState.rows }),
    });
    data = await res.json();
  } catch (e) {
    $("#primersSaveStatus").textContent = `저장 실패: ${e.message}`;
    toast("저장 실패: " + e.message, { type: "err", id: "primers-save" });
    return;
  }
  if (!res.ok) {
    primersState.errorsByRow = {};
    for (const e of data.errors || []) {
      primersState.errorsByRow[e.row] = primersState.errorsByRow[e.row] || {};
      primersState.errorsByRow[e.row][e.field] = e.message;
    }
    renderPrimersPage();
    $("#primersSaveStatus").textContent = `저장 실패: 입력값을 확인하세요 (${(data.errors || []).length}건)`;
    toast(`저장 실패: 입력값을 확인하세요 (${(data.errors || []).length}건)`, { type: "err", id: "primers-save" });
    return;
  }
  toast("프라이머 표를 저장했습니다.", { type: "ok", id: "primers-save" });
  primersState.rows = data.rows;
  primersState.savedIds = new Map(data.rows.filter((r) => r.primer_id).map((r) => [r.primer_id, r.primer_name]));
  if (!primersState.rows.length || primersState.rows[primersState.rows.length - 1].primer_name) {
    primersState.rows.push(primersBlankRow());
  }
  primersState.errorsByRow = {};
  renderPrimersPage();
  const purged = data.purged || [];
  $("#primersSaveStatus").textContent = purged.length
    ? `저장됨. 분석진행이 꺼진 ${purged.length}개 세트의 결과를 삭제했습니다 (${purged.map((x) => x.id).join(", ")}). 나머지는 "분석 다시하기"로 반영하세요.`
    : "저장됨. 반영하려면 \"분석 다시하기\"를 눌러주세요.";
  await Promise.all([loadSets(), loadGroups()]);                // a switched-off set disappears from the sidebar right away
  if (purged.length && !state.activeId) { /* overview/other pages keep their own refresh */ }
  if (state.activeType === "set" && purged.some((x) => x.id === state.activeId)) { state.activeId = null; state.activeType = null; renderSidebar(); }
}

$("#primersBtn").addEventListener("click", showPrimersPage);

// ---------------------------------------------------------------------------------------
// 표준물질 대표성 (output/standards/<seg>/representativeness.json, scripts/standard_representativeness.py)
// 실험실 양성대조물질(KDCA 표준물질)이 신뢰할 수 있는 전체 게놈을 얼마나 대표하는가 + 그 위치가
// 유전체 다양성 면에서 안정적인가. 등급: 합집합 일치율 >= 90% 적합, >= 50% 재검토 권고, 미만 심각.
const STD_GRADE = {
  good: { label: "적합", cls: "pill-ok" },
  review: { label: "재검토 권고", cls: "pill-warn" },
  critical: { label: "심각 · 표준물질 재선정 필요", cls: "pill-bad" },
  na: { label: "해당 없음", cls: "" },
};
const STD_LOCATION = {
  conserved: "보존 구간 (다양성 낮음)", moderate: "보통", variable: "고변이 구간 (다양성 높음)", unknown: "알 수 없음",
};
const STD_NA_REASON = {
  sites_outside_standard: "이 세트의 결합 부위가 표준물질 구간 밖에 있어 계산할 수 없습니다",
  no_binding_on_reference: "참조유전체에서 결합 부위를 찾지 못해 계산할 수 없습니다",
  unknown: "계산 결과가 없습니다",
};
let stdRepState = { data: null, seg: null };

function stdGradeBadge(grade, pct) {
  const g = STD_GRADE[grade] || STD_GRADE.na;
  return `<span class="pill ${g.cls}">${pct != null ? `${pct.toFixed(1)}% · ` : ""}${g.label}</span>`;
}

function stdAdvice(grade, union) {
  const weak = Object.entries(union.by_lineage || {})
    .filter(([, v]) => v.n >= 5 && v.matched / v.n < 0.5)
    .sort((a, b) => b[1].n - a[1].n)
    .map(([l, v]) => `${escHtml(l)} (${v.matched}/${v.n})`);
  const weakTxt = weak.length ? ` 대표하지 못하는 주요 리니지: ${weak.join(", ")}.` : "";
  if (grade === "good") return `표준물질이 전체 게놈의 대부분과 일치합니다. 그대로 사용해도 좋습니다.${weakTxt}`;
  if (grade === "review") return `일부 게놈과 일치하지 않습니다. 표준물질을 다시 고려하거나 보완 표준물질 추가를 검토하세요.${weakTxt}`;
  if (grade === "critical") return `절반 이상의 게놈이 표준물질과 일치하지 않습니다. 표준물질로 양성 확인이 되어도 실제 검출력을 보장하기 어려우니 재선정을 진지하게 고려하세요.${weakTxt}`;
  return "판정할 게놈이 없습니다.";
}

function renderStdRepSegment(seg, d, setNames) {
  const u = d.union;
  const p = d.params || {};
  const stdRows = Object.entries(d.standards).map(([id, info]) => {
    const ps = d.per_standard[id] || {};
    if (!info.mapped) return `<tr><td class="mono">${escHtml(id)}</td><td colspan="6">참조유전체에 정렬되지 않음</td></tr>`;
    const weakMap = info.mapped_fraction < 0.9;
    return `<tr>
      <td class="mono"><b>${escHtml(id)}</b></td>
      <td class="mono">${info.ref_start}–${info.ref_end}${info.reversed ? " <span class='note'>(역방향 보정)</span>" : ""}${weakMap ? ` <span class="pill pill-warn" title="표준물질 ${info.length}bp 중 ${(info.mapped_fraction * 100).toFixed(0)}%만 한 구간에 정렬됨">부분 정렬 ${(info.mapped_fraction * 100).toFixed(0)}%</span>` : ""}</td>
      <td>${STD_LOCATION[info.location] || info.location}<br><span class="note">평균 π ${info.mean_pi ?? "—"} · 백분위 ${info.percentile ?? "—"}</span></td>
      <td>${ps.matched ?? 0} / ${ps.assessable ?? 0}</td>
      <td>${ps.pct != null ? stdGradeBadge(ps.grade, ps.pct) : "—"}</td>
      <td>${ps.identity_median != null ? (ps.identity_median * 100).toFixed(1) + "%" : "—"}</td>
    </tr>`;
  }).join("");

  const lineageRows = Object.entries(u.by_lineage || {}).sort((a, b) => b[1].n - a[1].n).map(([l, v]) => {
    const pct = v.n ? (100 * v.matched) / v.n : 0;
    return `<div class="cov-row" style="grid-template-columns:110px 1fr 130px;cursor:default;">
      <div class="label">${escHtml(l)}${v.n < 10 ? ` <span class="note">n=${v.n}</span>` : ""}</div>
      <div class="cov-bar-track"><div class="cov-bar-fill validated" style="width:${pct}%"></div></div>
      <div class="cov-val" style="margin-left:0;">${v.matched} / ${v.n} (${pct.toFixed(0)}%)</div></div>`;
  }).join("");

  const sens = (u.sensitivity || []).map((s) =>
    `<span class="std-sens">identity ≥ ${(s.identity * 100).toFixed(0)}% → <b>${s.pct != null ? s.pct.toFixed(1) + "%" : "—"}</b></span>`).join(" ");

  const setRows = Object.entries(d.sets).map(([sid, e]) => {
    const nm = escHtml(setNames[sid] || sid);
    if (!e.applicable) {
      return `<tr class="genome-row" data-set-id="${sid}"><td>${nm}</td><td colspan="5" class="note">${STD_NA_REASON[e.reason] || STD_NA_REASON.unknown}</td></tr>`;
    }
    const locs = Object.entries(e.sites).map(([r, v]) => `${r}:${({ conserved: "보존", moderate: "보통", variable: "고변이" })[v.location] || "?"}`).join(" · ");
    const maj = Object.values(e.per_standard)[0]?.majority_agreement || {};
    const majTxt = Object.entries(maj).map(([r, v]) => `${r} ${(v * 100).toFixed(0)}%`).join(" · ");
    const sens2 = (e.sensitivity || []).map((s) => `${s.pct != null ? s.pct.toFixed(0) : "—"}%`).join(" / ");
    return `<tr class="genome-row" data-set-id="${sid}">
      <td>${nm}</td>
      <td class="mono">${e.covering_standards.map(escHtml).join(", ")}</td>
      <td>${e.union.matched} / ${e.union.assessable}</td>
      <td>${stdGradeBadge(e.grade, e.pct)}</td>
      <td>${locs}<br><span class="note">표준물질이 다수 염기와 같은 비율: ${majTxt}</span></td>
      <td class="note">${sens2}</td>
    </tr>`;
  }).join("");

  const unmatched = (u.unmatched_detail || []).map((r) =>
    `<tr><td class="mono">${escHtml(r.genome_id)}</td><td>${escHtml(r.lineage)}</td><td class="mono">${escHtml(r.best_standard || "—")}</td><td>${r.best_identity != null ? (r.best_identity * 100).toFixed(1) + "%" : "—"}</td></tr>`).join("");

  return `
    <div class="std-banner std-${u.grade}">
      <div class="std-banner-main">표준물질 전체(합집합) 대표성 ${stdGradeBadge(u.grade, u.pct)}</div>
      <div>${u.matched} / ${u.assessable}개 게놈이 최소 1개 표준물질과 일치${u.not_assessable ? ` (구간을 충분히 덮지 못한 ${u.not_assessable}개는 판정 제외)` : ""}</div>
      <div class="note" style="margin-top:4px;">${stdAdvice(u.grade, u)}</div>
    </div>
    <p class="note">"일치" 기준: 표준물질 구간의 ${(p.min_coverage * 100).toFixed(0)}% 이상을 덮는 게놈이 그 구간에서 <b>${(p.identity * 100).toFixed(0)}% 이상 동일</b>.
      등급 기준: ${p.good}% 이상 적합 · ${p.critical}% 이상 재검토 권고 · 미만 심각 (<code>config/config.yaml</code>의 <code>standard_rep_*</code>).
      기준 완화 시: ${sens}</p>

    <h3>표준물질별 결과</h3>
    <div class="table-wrap"><table class="genome-table"><thead><tr>
      <th>표준물질</th><th>참조유전체 위치</th><th>위치의 다양성 (pixy π)</th><th>일치 게놈</th><th>판정</th><th>일치 게놈의 중앙 identity</th>
    </tr></thead><tbody>${stdRows}</tbody></table></div>

    <h3>리니지별 일치 (전체 표준물질 합집합)</h3>
    <div class="cov-rows">${lineageRows}</div>

    <h3>프라이머 세트별: 결합 부위 기준 대표성</h3>
    <p class="note">이 세트의 F/R/Probe 결합 부위를 포함하는 표준물질만 대상으로, 부위의 염기가 표준물질과 같은지(F·R ≤ ${p.site_tol_primer}개, Probe ≤ ${p.site_tol_probe}개 차이) 비교합니다.
      마지막 열은 허용 차이를 (1/0)→(2/1)→(3/2)→(4/3)으로 늘렸을 때의 일치율입니다. 행을 클릭하면 세트 상세로 이동합니다.</p>
    <div class="table-wrap"><table class="genome-table"><thead><tr>
      <th>프라이머 세트</th><th>포함하는 표준물질</th><th>일치 게놈</th><th>판정</th><th>결합 부위의 다양성</th><th>허용 차이별 일치율</th>
    </tr></thead><tbody>${setRows}</tbody></table></div>

    <details style="margin-top:14px;"><summary>표준물질과 일치하지 않는 게놈 (상위 ${(u.unmatched_detail || []).length}개 표시)</summary>
      <p><a href="/api/standards/unmatched?seg=${seg}" download>⬇ 전체 목록 TSV 다운로드</a></p>
      <div class="table-wrap"><table class="genome-table"><thead><tr><th>게놈</th><th>리니지</th><th>가장 가까운 표준물질</th><th>identity</th></tr></thead><tbody>${unmatched}</tbody></table></div>
    </details>`;
}

function renderStdRepPage() {
  const view = $("#standardsView");
  const data = stdRepState.data;
  const segs = Object.keys(data.segments || {});
  if (!segs.length) {
    view.innerHTML = `<h2>표준물질 대표성</h2><p class="note">표준물질이 등록되어 있지 않거나 아직 계산되지 않았습니다. 표준물질은 선택 사항이며, "표준물질 입력 / 관리"에서 등록한 뒤 "분석 다시하기"를 실행하면 이 페이지가 채워집니다.</p>`;
    return;
  }
  if (!stdRepState.seg || !segs.includes(stdRepState.seg)) stdRepState.seg = segs[0];
  const seg = stdRepState.seg;
  view.innerHTML = `
    <h2>표준물질 대표성 (KDCA 양성대조물질)</h2>
    <p class="note">표준물질이 신뢰할 수 있는 전체 LASV 게놈(품질 필터 통과 게놈)을 얼마나 대표하는지, 그리고 표준물질/프라이머 결합 부위가 유전체 다양성 면에서 안정적인 곳에 있는지 확인합니다.
      표준물질에서 양성이 나와도 실제 바이러스 대부분이 표준물질과 크게 다르면 그 양성은 실제 검출력을 보장하지 못합니다.</p>
    <div class="seg-tab-row">${segs.map((s) => `<button class="seg-tab ${s === seg ? "active" : ""}" data-seg="${s}">${s.toUpperCase()} segment (${data.segments[s].genome_count}개 게놈)</button>`).join("")}</div>
    <div id="stdRepBody">${renderStdRepSegment(seg, data.segments[seg], data.set_names || {})}</div>`;
  view.querySelectorAll(".seg-tab").forEach((b) => b.addEventListener("click", () => { stdRepState.seg = b.dataset.seg; renderStdRepPage(); }));
  view.querySelectorAll("tr[data-set-id]").forEach((tr) => tr.addEventListener("click", () => selectSet(tr.dataset.setId)));
}

async function showStandardsPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("standardsView");
  const view = $("#standardsView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    stdRepState.data = await fetchJSON("/api/standards/representativeness");
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
    return;
  }
  renderStdRepPage();
}

$("#standardsBtn").addEventListener("click", showStandardsPage);

// ---------------------------------------------------------------------------------------
// 메타데이터 표: every analysed population genome (S + L) as a spreadsheet -- Nextstrain metadata,
// how many primer sets validate it, diversity-QC verdict. Click a column title to sort
// (asc -> desc -> off); every column has a Numbers-style Quick Filter (funnel: tick the values to
// keep) plus a typed filter (text: contains / "a|b" = either / "=x" exact; numbers: ">=100",
// "<5", "100-200"); a global search box and S/L toggle sit on top.
const mtState = { cols: null, rows: null, sort: { key: null, dir: 1 }, filters: {}, global: "", seg: "ALL", view: [], menu: null, passedOnly: false };
const MT_CAP = 5000;
const MT_PASSED = "통과";   // data_qc value of a genome that passed the length/N filter

function mtFmt(col, v) {
  if (v === null || v === undefined || v === "") return "";
  if (col.key === "identity_to_ref" || col.key === "coverage") return Number(v).toFixed(3);
  if (col.key === "validated_pct") return Number(v).toFixed(1);
  return String(v);
}

function mtNumPredicate(expr) {
  const e = expr.trim();
  if (!e) return null;
  const range = e.match(/^(-?\d+(?:\.\d+)?)\s*(?:-|~|\.\.)\s*(-?\d+(?:\.\d+)?)$/);
  if (range) { const lo = +range[1], hi = +range[2]; return (v) => v !== null && v >= lo && v <= hi; }
  const tokens = e.split(/[,\s]+/).filter(Boolean);
  const preds = tokens.map((t) => {
    const m = t.match(/^(>=|<=|>|<|=)?(-?\d+(?:\.\d+)?)$/);
    if (!m) return null;
    const n = +m[2], op = m[1] || "=";
    return (v) => v !== null && (op === ">=" ? v >= n : op === "<=" ? v <= n : op === ">" ? v > n : op === "<" ? v < n : v === n);
  });
  if (preds.every(Boolean)) return (v) => preds.every((p) => p(v));
  const needle = e.toLowerCase();
  return (v) => v !== null && String(v).toLowerCase().includes(needle);
}

function mtTextPredicate(expr) {
  const e = expr.trim();
  if (!e) return null;
  if (e.startsWith("=")) { const x = e.slice(1).trim().toLowerCase(); return (v) => String(v ?? "").toLowerCase() === x; }
  const alts = e.split("|").map((x) => x.trim().toLowerCase()).filter(Boolean);
  return (v) => { const s = String(v ?? "").toLowerCase(); return alts.some((a) => s.includes(a)); };
}

function mtRowPasses(row, skipKey) {
  const g = mtState.global.trim().toLowerCase();
  if (mtState.seg !== "ALL" && row.segment !== mtState.seg) return false;
  if (mtState.passedOnly && row.data_qc !== MT_PASSED) return false;
  if (g && !mtState.cols.some((c) => String(row[c.key] ?? "").toLowerCase().includes(g))) return false;
  for (const c of mtState.cols) {
    if (c.key === skipKey) continue;
    const f = mtState.filters[c.key];
    if (!f) continue;
    const v = row[c.key];
    if (f.set && !f.set.has(v === null || v === undefined || v === "" ? "" : String(v))) return false;
    if (f.pred && !f.pred(v === undefined ? null : v)) return false;
  }
  return true;
}

function mtCompute() {
  const rows = mtState.rows.filter((r) => mtRowPasses(r, null));
  const { key, dir } = mtState.sort;
  if (key) {
    const col = mtState.cols.find((c) => c.key === key);
    rows.sort((a, b) => {
      const av = a[key], bv = b[key];
      const an = av === null || av === undefined || av === "", bn = bv === null || bv === undefined || bv === "";
      if (an || bn) return an && bn ? 0 : an ? 1 : -1;   // blanks always last
      return dir * (col.type === "num" ? av - bv : String(av).localeCompare(String(bv), undefined, { numeric: true }));
    });
  }
  mtState.view = rows;
}

function mtRenderBody() {
  mtCompute();
  const { cols, view } = mtState;
  const tbody = $("#mtBody");
  const shown = view.slice(0, MT_CAP);
  tbody.innerHTML = shown.map((r, i) => `<tr data-i="${i}" class="${r.data_qc === MT_PASSED ? "mt-passed" : "mt-dim"}">${cols.map((c, ci) => {
    const v = mtFmt(c, r[c.key]);
    const cls = [c.type === "num" ? "mt-num" : "", ci === 0 ? "mt-sticky mt-s0" : ci === 1 ? "mt-sticky mt-s1" : ""].join(" ").trim();
    let inner = escHtml(v);
    if (c.key === "segment") inner = `<span class="pill pill-${v}">${escHtml(v)}</span>`;
    else if (c.key === "validated_pct" && v !== "") inner = `<span class="mt-bar" style="width:${Math.min(100, r[c.key])}%"></span><span class="mt-bar-text">${escHtml(v)}</span>`;
    else if (c.key === "data_qc") inner = `<span class="mt-qc ${v === MT_PASSED ? "ok" : "bad"}">${v === MT_PASSED ? "✓ " : ""}${escHtml(v)}</span>`;
    else if (c.key === "diversity_qc") inner = `<span class="mt-qc ${v === "사용" ? "ok" : v.startsWith("제외") ? "bad" : ""}">${escHtml(v)}</span>`;
    return `<td class="${cls}" title="${escHtml(v)}">${inner}</td>`;
  }).join("")}</tr>`).join("");
  $("#mtCount").textContent = `${view.length.toLocaleString()} / ${mtState.rows.length.toLocaleString()}행` + (view.length > MT_CAP ? ` (처음 ${MT_CAP.toLocaleString()}행만 표시)` : "");
  // header state: sort arrows + active-filter highlight
  cols.forEach((c) => {
    const th = document.querySelector(`#mtHead th[data-key="${c.key}"]`);
    if (!th) return;
    th.classList.toggle("sorted", mtState.sort.key === c.key);
    th.querySelector(".mt-arrow").textContent = mtState.sort.key === c.key ? (mtState.sort.dir === 1 ? "▲" : "▼") : "";
    const f = mtState.filters[c.key];
    th.querySelector(".mt-funnel").classList.toggle("on", !!(f && f.set));
  });
  mtRenderChips();
}

function mtRenderChips() {
  const chips = [];
  for (const c of mtState.cols) {
    const f = mtState.filters[c.key];
    if (f && f.set) chips.push(`<span class="mt-chip" data-clear="${c.key}">${escHtml(c.label)}: ${f.set.size}개 값 <b>×</b></span>`);
  }
  $("#mtChips").innerHTML = chips.join("");
}

function mtCloseMenu() {
  if (mtState.menu) { mtState.menu.remove(); mtState.menu = null; }
  document.removeEventListener("mousedown", mtOutside, true);
  document.removeEventListener("keydown", mtMenuKey, true);
}
function mtOutside(ev) { if (mtState.menu && !mtState.menu.contains(ev.target) && !ev.target.closest(".mt-funnel")) mtCloseMenu(); }
function mtMenuKey(ev) { if (ev.key === "Escape") { ev.stopPropagation(); mtCloseMenu(); } }

// Numbers-style Quick Filter: tick the values to keep. Counts reflect the other active filters.
function mtOpenQuickFilter(key, anchor) {
  mtCloseMenu();
  const col = mtState.cols.find((c) => c.key === key);
  const counts = new Map();
  for (const r of mtState.rows) {
    if (!mtRowPasses(r, key)) continue;
    const v = r[key] === null || r[key] === undefined ? "" : String(r[key]);
    counts.set(v, (counts.get(v) || 0) + 1);
  }
  const all = [...counts.keys()].sort((a, b) => col.type === "num" && a !== "" && b !== "" ? a - b : (a === "" ? 1 : b === "" ? -1 : a.localeCompare(b, undefined, { numeric: true })));
  const cur = mtState.filters[key]?.set;
  const selected = new Set(cur ? [...cur] : all);
  const menu = el("div", "mt-menu");
  const rect = anchor.getBoundingClientRect();
  menu.style.left = `${Math.max(8, Math.min(rect.left, window.innerWidth - 270))}px`;
  menu.style.top = `${rect.bottom + 4}px`;
  menu.innerHTML = `<div class="mt-menu-head">${escHtml(col.label)} 필터</div>
    <input type="text" class="mt-menu-search" placeholder="값 검색..." spellcheck="false">
    <div class="mt-menu-actions"><a href="#" data-act="all">모두 선택</a><a href="#" data-act="none">모두 해제</a></div>
    <div class="mt-menu-list"></div>`;
  document.body.appendChild(menu);
  mtState.menu = menu;
  const listEl = menu.querySelector(".mt-menu-list");
  const draw = () => {
    const q = menu.querySelector(".mt-menu-search").value.trim().toLowerCase();
    const vis = all.filter((v) => !q || (v === "" ? "(비어 있음)" : v).toLowerCase().includes(q)).slice(0, 300);
    listEl.innerHTML = vis.map((v) => `<label class="mt-menu-item"><input type="checkbox" data-v="${escHtml(v)}" ${selected.has(v) ? "checked" : ""}>
      <span class="mt-menu-val">${v === "" ? "<i>(비어 있음)</i>" : escHtml(v)}</span><span class="mt-menu-n">${counts.get(v)}</span></label>`).join("")
      || "<div class='note' style='padding:6px 8px;'>일치하는 값 없음</div>";
  };
  const apply = () => {
    if (selected.size === all.length) delete mtState.filters[key]?.set;
    else (mtState.filters[key] = mtState.filters[key] || {}).set = new Set(selected);
    if (mtState.filters[key] && !mtState.filters[key].set && !mtState.filters[key].pred) delete mtState.filters[key];
    mtRenderBody();
  };
  listEl.addEventListener("change", (ev) => {
    const cb = ev.target.closest("input[data-v]");
    if (!cb) return;
    const v = cb.dataset.v;
    if (cb.checked) selected.add(v); else selected.delete(v);
    apply();
  });
  menu.querySelector(".mt-menu-search").addEventListener("input", draw);
  menu.querySelectorAll(".mt-menu-actions a").forEach((a) => a.addEventListener("click", (ev) => {
    ev.preventDefault();
    const q = menu.querySelector(".mt-menu-search").value.trim().toLowerCase();
    const targets = all.filter((v) => !q || (v === "" ? "(비어 있음)" : v).toLowerCase().includes(q));
    targets.forEach((v) => (a.dataset.act === "all" ? selected.add(v) : selected.delete(v)));
    draw(); apply();
  }));
  draw();
  menu.querySelector(".mt-menu-search").focus();
  document.addEventListener("mousedown", mtOutside, true);
  document.addEventListener("keydown", mtMenuKey, true);
}

function mtExportCsv() {
  const { cols, view } = mtState;
  const q = (v) => { const s = v === null || v === undefined ? "" : String(v); return /[",\n\t]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
  const text = [cols.map((c) => q(c.label)).join(","), ...view.map((r) => cols.map((c) => q(r[c.key])).join(","))].join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob(["﻿" + text], { type: "text/csv;charset=utf-8" }));
  a.download = "lasv_metadata_table.csv";
  a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
}

function mtResetFilters() {
  mtState.filters = {}; mtState.global = ""; mtState.sort = { key: null, dir: 1 }; mtState.passedOnly = false;
  $("#mtGlobal").value = ""; $("#mtPassedOnly").checked = false;
  document.querySelectorAll("#mtFilterRow input").forEach((i) => { i.value = ""; });
  mtRenderBody();
}

function renderMetaTablePage() {
  const view = $("#metaTableView");
  const { cols } = mtState;
  const head = cols.map((c, i) => `<th data-key="${c.key}" class="${c.type === "num" ? "mt-num" : ""} ${i === 0 ? "mt-sticky mt-s0" : i === 1 ? "mt-sticky mt-s1" : ""}">
      <div class="mt-th"><span class="mt-title" title="클릭: 정렬 (오름차순 → 내림차순 → 해제)">${escHtml(c.label)}</span><span class="mt-arrow"></span>
      <button type="button" class="mt-funnel" title="퀵 필터: 값 선택">&#9662;</button></div></th>`).join("");
  const filt = cols.map((c, i) => `<th class="${i === 0 ? "mt-sticky mt-s0" : i === 1 ? "mt-sticky mt-s1" : ""}"><input type="text" data-key="${c.key}" placeholder="${c.type === "num" ? "예: >=100" : "포함…"}" spellcheck="false" autocomplete="off"></th>`).join("");
  view.innerHTML = `
    <h2>메타데이터 표</h2>
    <p class="note" style="margin-top:4px;">분석에 쓰인 population 게놈(S+L)의 Nextstrain 메타데이터와 검증 결과를 표로 봅니다. 열 제목을 누르면 정렬되고, 제목 옆 ▾는 숫자/값을 골라 거르는
      퀵 필터, 아래 입력칸은 직접 입력 필터입니다 (텍스트: 포함 · <code>a|b</code> 둘 중 하나 · <code>=정확히</code> / 숫자: <code>&gt;=1000</code> · <code>&lt;5</code> · <code>100-200</code>).</p>
    <div class="mt-toolbar">
      <input id="mtGlobal" type="search" class="mt-global" placeholder="전체 검색 (모든 열)…" spellcheck="false">
      <div class="mt-seg">${["ALL", "S", "L"].map((s) => `<button type="button" data-seg="${s}" class="${s === mtState.seg ? "active" : ""}">${s === "ALL" ? "전체" : s}</button>`).join("")}</div>
      <label class="mt-check" title="길이·모호염기(N)·중복 필터를 통과해 분석에 쓰인 게놈만 표시"><input type="checkbox" id="mtPassedOnly" ${mtState.passedOnly ? "checked" : ""}> 데이터 QC 통과만</label>
      <button id="mtReset" class="btn btn-secondary">필터·정렬 초기화</button>
      <button id="mtCsv" class="btn btn-secondary" title="현재 보이는 행(필터/정렬 적용)을 CSV로 저장">&#8681; CSV</button>
      <span id="mtCount" class="note mt-count"></span>
    </div>
    <div id="mtSummary" class="note mt-summary"></div>
    <div id="mtChips" class="mt-chips"></div>
    <div class="mt-wrap"><table class="mt-table">
      <thead id="mtHead"><tr>${head}</tr><tr id="mtFilterRow" class="mt-filter-row">${filt}</tr></thead>
      <tbody id="mtBody"></tbody></table></div>`;

  $("#mtGlobal").addEventListener("input", (ev) => { mtState.global = ev.target.value; mtDebounced(); });
  view.querySelectorAll(".mt-seg button").forEach((b) => b.addEventListener("click", () => {
    mtState.seg = b.dataset.seg;
    view.querySelectorAll(".mt-seg button").forEach((x) => x.classList.toggle("active", x === b));
    mtRenderBody();
  }));
  $("#mtPassedOnly").addEventListener("change", (ev) => { mtState.passedOnly = ev.target.checked; mtRenderBody(); });
  const passed = mtState.rows.filter((r) => r.data_qc === MT_PASSED);
  const bySeg = (sg) => passed.filter((r) => r.segment === sg).length;
  const div = passed.filter((r) => r.diversity_qc === "사용").length;
  $("#mtSummary").innerHTML = `전체 <b>${mtState.rows.length.toLocaleString()}</b>개 중 <b class="mt-qc ok">데이터 QC 통과 ${passed.length.toLocaleString()}개</b> (S ${bySeg("S")} · L ${bySeg("L")}) — 길이·모호염기(N)·중복 필터를 통과해 분석에 쓰인 게놈입니다. 나머지(흐리게 표시)는 제외된 사유가 "데이터 QC" 열에 나옵니다. 그중 diversity 품질 필터까지 통과해 π 계산에 쓰인 게놈은 ${div.toLocaleString()}개입니다.`;
  $("#mtReset").addEventListener("click", mtResetFilters);
  $("#mtCsv").addEventListener("click", mtExportCsv);
  $("#mtHead").addEventListener("click", (ev) => {
    const funnel = ev.target.closest(".mt-funnel");
    const th = ev.target.closest("th[data-key]");
    if (funnel && th) { mtOpenQuickFilter(th.dataset.key, funnel); return; }
    if (th && ev.target.closest(".mt-title, .mt-arrow")) {
      const key = th.dataset.key;
      const s = mtState.sort;
      if (s.key !== key) mtState.sort = { key, dir: 1 };
      else if (s.dir === 1) mtState.sort = { key, dir: -1 };
      else mtState.sort = { key: null, dir: 1 };
      mtRenderBody();
    }
  });
  $("#mtFilterRow").addEventListener("input", (ev) => {
    const inp = ev.target.closest("input[data-key]");
    if (!inp) return;
    const key = inp.dataset.key;
    const col = cols.find((c) => c.key === key);
    const pred = col.type === "num" ? mtNumPredicate(inp.value) : mtTextPredicate(inp.value);
    const f = (mtState.filters[key] = mtState.filters[key] || {});
    f.pred = pred;
    if (!f.pred && !f.set) delete mtState.filters[key];
    mtDebounced();
  });
  $("#mtChips").addEventListener("click", (ev) => {
    const chip = ev.target.closest("[data-clear]");
    if (!chip) return;
    const f = mtState.filters[chip.dataset.clear];
    if (f) { delete f.set; if (!f.pred) delete mtState.filters[chip.dataset.clear]; }
    mtRenderBody();
  });
  mtRenderBody();
}

let _mtTimer = null;
function mtDebounced() { clearTimeout(_mtTimer); _mtTimer = setTimeout(mtRenderBody, 120); }

async function showMetaTablePage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("metaTableView");
  const view = $("#metaTableView");
  if (!mtState.rows) {
    view.innerHTML = "<p class='note'>불러오는 중...</p>";
    try {
      const data = await fetchJSON("/api/metadata/table");
      mtState.cols = data.columns; mtState.rows = data.rows;
    } catch (e) {
      view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
      return;
    }
  }
  renderMetaTablePage();
}

$("#metaTableBtn").addEventListener("click", showMetaTablePage);

// ---------------------------------------------------------------------------------------
// 요구사항 기록: everything asked for since the project started, with what was reviewed/changed and the
// status. Read-only view of docs/requirements_log.json (served by /api/requirements): summary counts,
// text search, category and status filters, one table per phase.
const REQ_STATUS = { done: "완료", partial: "일부 미검증", answer: "답변·제안", dropped: "철회·불필요", open: "진행 중" };
let reqState = { data: null, q: "", cat: "", st: "" };

async function showRequirementsPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("requirementsView");
  const view = $("#requirementsView");
  if (!reqState.data) {
    view.innerHTML = "<p class='note'>불러오는 중...</p>";
    try { reqState.data = await fetchJSON("/api/requirements"); }
    catch (e) { view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`; return; }
  }
  renderRequirements();
}

function renderRequirements() {
  const d = reqState.data, view = $("#requirementsView");
  const all = d.phases.flatMap((p) => p.items);
  if (!all.length) { view.innerHTML = "<h2>요구사항 기록</h2><p class='note'>기록이 없습니다.</p>"; return; }
  const counts = {}; all.forEach((i) => { counts[i.status] = (counts[i.status] || 0) + 1; });
  const cats = [...new Set(all.map((i) => i.cat))].sort();
  const code = (t) => escHtml(t || "").replace(/`([^`]+)`/g, "<code>$1</code>");
  let n = 0;
  const body = d.phases.map((p) => {
    const rows = p.items.map((i) => {
      n += 1;
      return `<tr data-cat="${escHtml(i.cat)}" data-st="${i.status}"><td class="rq-n">${n}</td><td class="rq-d">${escHtml(i.date)}</td><td><span class="rq-cat">${escHtml(i.cat)}</span></td>
        <td class="rq-req">${code(i.req)}</td><td>${code(i.rev)}${i.note ? `<div class="rq-note">${code(i.note)}</div>` : ""}</td><td><span class="rq-st rq-${i.status}">${REQ_STATUS[i.status] || i.status}</span></td></tr>`;
    }).join("");
    return `<tr class="rq-phase"><td colspan="6">${escHtml(p.title)}</td></tr>${rows}`;
  }).join("");
  const sumTiles = Object.keys(REQ_STATUS).filter((k) => counts[k]).map((k) => `<div class="rq-tile"><b class="rq-${k}-t">${counts[k]}</b><span>${REQ_STATUS[k]}</span></div>`).join("");
  view.innerHTML = `
    <h2>요구사항 기록</h2>
    <p class="note" style="margin-top:2px;">프로젝트 시작(2026-09-14)부터 요청하신 사항과 검토·수정 내용, 처리 상태입니다. 갱신일 ${escHtml(d.updated || "")} · 데이터는 <code>docs/requirements_log.json</code>.</p>
    <div class="rq-sum"><div class="rq-tile"><b>${all.length}</b><span>전체 요구</span></div>${sumTiles}</div>
    <div class="rq-tools"><input id="rqSearch" type="search" placeholder="요구 사항·수정 내용 검색" spellcheck="false" value="${escHtml(reqState.q)}">
      <span class="rq-chips" id="rqCats"><button type="button" data-cat="" class="rq-chip ${reqState.cat === "" ? "on" : ""}">전체</button>${cats.map((c) => `<button type="button" data-cat="${escHtml(c)}" class="rq-chip ${reqState.cat === c ? "on" : ""}">${escHtml(c)}</button>`).join("")}</span>
      <select id="rqStatus"><option value="">모든 상태</option>${Object.entries(REQ_STATUS).filter(([k]) => counts[k]).map(([k, v]) => `<option value="${k}" ${reqState.st === k ? "selected" : ""}>${v}</option>`).join("")}</select></div>
    <div class="rq-wrap"><table class="rq-table"><thead><tr><th>#</th><th>날짜</th><th>구분</th><th>요구 사항</th><th>검토 · 수정 내용</th><th>상태</th></tr></thead><tbody id="rqBody">${body}</tbody></table></div>
    <h3 style="margin-top:22px;">남은 확인 사항</h3>
    <ul class="rq-open">${(d.open || []).map((o) => `<li><b>${code(o.item)}</b><span>${code(o.note)}</span></li>`).join("")}</ul>`;
  const apply = () => {
    const q = reqState.q.trim().toLowerCase(); let vis = false;
    const trs = [...document.querySelectorAll("#rqBody tr")];
    for (let i = trs.length - 1; i >= 0; i--) {
      const r = trs[i];
      if (r.classList.contains("rq-phase")) { r.hidden = !vis; vis = false; continue; }
      const ok = (!reqState.cat || r.dataset.cat === reqState.cat) && (!reqState.st || r.dataset.st === reqState.st) && (!q || r.textContent.toLowerCase().includes(q));
      r.hidden = !ok; if (ok) vis = true;
    }
  };
  $("#rqSearch").addEventListener("input", (ev) => { reqState.q = ev.target.value; apply(); });
  $("#rqStatus").addEventListener("change", (ev) => { reqState.st = ev.target.value; apply(); });
  $("#rqCats").addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-cat]"); if (!b) return;
    reqState.cat = b.dataset.cat;
    document.querySelectorAll("#rqCats .rq-chip").forEach((x) => x.classList.toggle("on", x === b));
    apply();
  });
  apply();
}

$("#requirementsBtn").addEventListener("click", showRequirementsPage);

// ---------------------------------------------------------------------------------------
// 표준물질 입력 / 관리: the reference-standard (positive-control) sequences as a spreadsheet, like the
// primer page. Saved to config/standards.tsv; the server regenerates the FASTA/metadata files the
// pipeline reads. The NAME is the standard's id in every result and wet-lab record, so renaming or
// deleting a saved one asks for confirmation. Long sequences are edited in a window (click the cell),
// and several can be added at once by pasting FASTA.
let stdMgr = { rows: [], saved: new Map(), errors: {}, general: [] };

function stdBlankRow() { return { name: "", segment: "S", sequence: "", source: "", remark: "", proceed: "True", add_date: "" }; }
const stdClean = (s) => (s || "").replace(/\s+/g, "").toUpperCase();

async function showStdManagerPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("stdMgrView");
  const view = $("#stdMgrView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    const data = await fetchJSON("/api/standards/manifest");
    stdMgr.rows = data.rows || [];
    stdMgr.saved = new Map(stdMgr.rows.map((r) => [r.name, r.name]));
    stdMgr.rows.push(stdBlankRow());
    stdMgr.errors = {}; stdMgr.general = [];
    renderStdManager();
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`;
  }
}

function stdRowHtml(r, i) {
  const errs = stdMgr.errors[i] || {};
  const ec = (f) => (errs[f] ? "cell-error" : "");
  const txt = (f, ph = "") => `<td class="${ec(f)}"><input type="text" class="cell-input" data-field="${f}" value="${escHtml(r[f] || "")}" placeholder="${ph}" spellcheck="false" autocomplete="off"></td>`;
  const seq = stdClean(r.sequence);
  const on = ["true", "1", "yes", "y"].includes(String(r.proceed).toLowerCase());
  const errRow = Object.keys(errs).length
    ? `<tr class="err-row"><td class="row-num"></td><td colspan="9"><div class="primers-cell-errors">${Object.values(errs).map(escHtml).join(" / ")}</div></td></tr>` : "";
  return `<tr data-row="${i}" class="${r.add_date ? "" : "row-new"}">
    <td class="row-num">${i + 1}</td>
    ${txt("name", "예: NP1")}
    <td class="${ec("segment")}"><select class="cell-input" data-field="segment"><option value="S" ${r.segment === "S" ? "selected" : ""}>S</option><option value="L" ${r.segment === "L" ? "selected" : ""}>L</option></select></td>
    <td class="std-seq-cell ${ec("sequence")}"><button type="button" class="std-seq-btn" data-edit-seq="${i}" title="클릭해서 서열 편집·붙여넣기">${seq ? `<span class="mono">${escHtml(seq.slice(0, 22))}${seq.length > 22 ? "…" : ""}</span> <span class="std-len">${seq.length.toLocaleString()} bp</span>` : "<i>서열 입력…</i>"}</button></td>
    ${txt("source", "예: KDCA 2024")}
    ${txt("remark", "")}
    <td class="center-cell"><input type="checkbox" class="cell-input" data-field="proceed" ${on ? "checked" : ""} title="해제하면 분석에서 제외됩니다"></td>
    <td class="add-date-cell">${escHtml(r.add_date || "저장 시 자동")}</td>
    <td class="del-cell"><button class="row-del-btn" data-del-row="${i}" title="행 삭제">&times;</button></td>
  </tr>${errRow}`;
}

function stdXlColumns() {
  return [
    { key: "name", title: "이름 (ID)", width: 190, sticky: true },
    { key: "segment", title: "Seg", width: 52, upper: true },
    { key: "sequence", title: "서열", width: 420, seq: true },
    { key: "length", title: "길이 (bp)", width: 80, ro: true, num: true, value: (r) => stdClean(r.sequence).length || "" },
    { key: "source", title: "출처", width: 150 },
    { key: "remark", title: "비고", width: 190 },
    { key: "proceed", title: "분석진행", width: 84, type: "bool" },
    { key: "add_date", title: "추가일", width: 100, ro: true, value: (r) => r.add_date || "저장 시 자동" },
  ];
}

function renderStdManager() {
  const view = $("#stdMgrView");
  const xlOn = xlModeOn("standards");
  const gen = stdMgr.general.length ? `<div class="primers-cell-errors" style="margin:8px 0;">${stdMgr.general.map(escHtml).join("<br>")}</div>` : "";
  view.innerHTML = `
    <h2>표준물질 입력 / 관리</h2>
    <div class="note">
      실험실 양성대조물질(표준물질)의 서열을 등록합니다. 저장하면 분석에 쓰이는 파일이 자동으로 만들어지고, 반영하려면 사이드바의 "분석 다시하기"를 눌러야 합니다.
      <b>이름은 결과·실험 검증 기록의 ID</b>로 쓰이므로 저장된 이름을 바꾸거나 지우면 기존 결과·기록과 연결이 끊깁니다. <b>분석진행</b>을 해제한 행은 분석에서 제외됩니다.
      표준물질은 선택 사항입니다 — 등록하지 않으면 표준물질 평가는 건너뜁니다.
    </div>
    <div class="primers-toolbar">
      <button id="stdAddRowBtn" class="btn btn-secondary">+ 행 추가</button>
      <button id="stdXlBtn" class="btn btn-secondary xl-mode-btn${xlOn ? " is-on" : ""}" title="엑셀처럼 셀을 선택·복사·붙여넣기·실행 취소할 수 있는 격자 화면으로 전환">스프레드시트 형식</button>
      <button id="stdFastaBtn" class="btn btn-secondary" title="FASTA 텍스트를 붙여넣어 여러 개를 한 번에 추가">FASTA로 추가</button>
      <a id="stdXlsxDl" class="btn btn-secondary" href="/api/standards/xlsx" download title="현재 표준물질 목록을 엑셀(.xlsx)로 내려받습니다 (열 설명 시트 포함)">&#11015; 엑셀 다운로드</a>
      <button id="stdXlsxUp" class="btn btn-secondary" title="엑셀(.xlsx) 파일의 행을 이 표로 불러옵니다 (저장 전에 확인할 수 있습니다)">&#11014; 엑셀 업로드</button>
      <input type="file" id="stdXlsxFile" accept=".xlsx" hidden>
      <button id="stdSaveBtn" class="btn btn-primary">저장</button>
      <span id="stdSaveStatus" class="form-status"></span>
    </div>${gen}
    ${xlOn ? '<div id="stdXl" class="xl-host"></div>' : ""}
    <div class="primers-grid-wrap" ${xlOn ? "hidden" : ""}><table class="primers-table std-table">
      <thead><tr><th class="row-num">#</th><th>이름 (ID)</th><th>Seg</th><th>서열</th><th>출처</th><th>비고</th><th style="text-align:center;">분석진행</th><th>추가일</th><th></th></tr></thead>
      <tbody id="stdBody">${stdMgr.rows.map(stdRowHtml).join("")}</tbody></table></div>
    <div class="primers-hint">서열은 IUPAC 코드만 허용되며 공백·줄바꿈은 자동으로 제거됩니다. 서열 칸을 누르면 큰 편집 창이 열립니다.</div>`;
  $("#stdXlBtn").addEventListener("click", () => { xlModeSet("standards", !xlOn); renderStdManager(); });
  $("#stdFastaBtn").addEventListener("click", stdImportFasta);
  wireXlsxButtons("standards", $("#stdXlsxUp"), $("#stdXlsxFile"), {
    existing: () => stdMgr.rows.filter((r) => r.name || r.sequence),
    columns: stdXlColumns, blankRow: stdBlankRow,
    apply: (rows, mode) => {
      const base = mode === "replace" ? [] : stdMgr.rows.filter((r) => r.name || r.sequence);
      const add = rows.map((r) => (mode === "append" ? { ...r, add_date: "" } : r));
      stdMgr.rows = [...base, ...add, stdBlankRow()];
      stdMgr.errors = {}; stdMgr.general = [];
      renderStdManager();
      $("#stdSaveStatus").textContent = `엑셀에서 ${rows.length}행을 불러왔습니다 — 확인 후 저장을 눌러 반영하세요.`;
    },
  });
  $("#stdSaveBtn").addEventListener("click", stdSave);
  if (xlOn) {
    $("#stdAddRowBtn").addEventListener("click", () => { stdMgr.rows.push(stdBlankRow()); renderStdManager(); });
    createXlGrid({
      host: $("#stdXl"), columns: stdXlColumns(), getRows: () => stdMgr.rows, blankRow: stdBlankRow,
      getErrors: () => stdMgr.errors,
      ensureBlank: () => { const l = stdMgr.rows[stdMgr.rows.length - 1]; if (!l || l.name || l.sequence) stdMgr.rows.push(stdBlankRow()); },
      onChange: () => { stdMgr.errors = {}; },
    });
    return;
  }
  const body = $("#stdBody");
  const ensureBlank = () => {
    const last = stdMgr.rows[stdMgr.rows.length - 1];
    if (!last || last.name || last.sequence) { stdMgr.rows.push(stdBlankRow()); return true; }
    return false;
  };
  body.addEventListener("input", (ev) => {
    const inp = ev.target.closest("[data-field]"); const tr = ev.target.closest("tr[data-row]");
    if (!inp || !tr) return;
    const i = +tr.dataset.row;
    stdMgr.rows[i][inp.dataset.field] = inp.type === "checkbox" ? (inp.checked ? "True" : "False") : inp.value;
    if (i === stdMgr.rows.length - 1 && inp.dataset.field === "name" && inp.value && ensureBlank()) body.insertAdjacentHTML("beforeend", stdRowHtml(stdMgr.rows[stdMgr.rows.length - 1], stdMgr.rows.length - 1));
  });
  body.addEventListener("click", (ev) => {
    const ed = ev.target.closest("[data-edit-seq]");
    if (ed) return stdEditSequence(+ed.dataset.editSeq);
    const del = ev.target.closest("[data-del-row]");
    if (del) { stdMgr.rows.splice(+del.dataset.delRow, 1); ensureBlank(); stdMgr.errors = {}; renderStdManager(); }
  });
  $("#stdAddRowBtn").addEventListener("click", () => { stdMgr.rows.push(stdBlankRow()); renderStdManager(); });
}

function stdEditSequence(i) {
  const r = stdMgr.rows[i];
  const m = openDetailModal({
    title: `서열 편집 <span class="mono">${escHtml(r.name || "(이름 없음)")}</span>`, subtitle: `segment ${escHtml(r.segment)}`,
    html: `<textarea id="stdSeqText" class="std-seq-text mono" spellcheck="false" placeholder="서열을 붙여넣으세요 (FASTA 한 개도 가능)">${escHtml(r.sequence || "")}</textarea>
      <div class="ui-dialog-actions" style="padding:10px 0 0;"><span id="stdSeqInfo" class="note" style="margin:0 auto 0 0;"></span>
      <button type="button" class="btn btn-secondary" id="stdSeqCancel">취소</button><button type="button" class="btn btn-primary" id="stdSeqOk">적용</button></div>`,
  });
  const ta = m.overlay.querySelector("#stdSeqText"), info = m.overlay.querySelector("#stdSeqInfo");
  const parse = () => ta.value.split("\n").filter((l) => !l.startsWith(">")).join("");
  const upd = () => { const s = stdClean(parse()); const bad = [...new Set(s.replace(/[ACGTURYSWKMBDHVN]/g, ""))]; info.textContent = `${s.length.toLocaleString()} bp` + (bad.length ? ` · 올바르지 않은 문자: ${bad.join(" ")}` : ""); info.style.color = bad.length ? "var(--bad)" : ""; };
  ta.addEventListener("input", upd); upd(); ta.focus();
  m.overlay.querySelector("#stdSeqCancel").addEventListener("click", m.close);
  m.overlay.querySelector("#stdSeqOk").addEventListener("click", () => { r.sequence = stdClean(parse()); stdMgr.errors = {}; m.close(); renderStdManager(); });
}

function stdImportFasta() {
  const m = openDetailModal({
    title: "FASTA로 표준물질 추가", subtitle: "레코드 하나가 표준물질 하나가 됩니다 (헤더의 첫 단어가 이름)",
    html: `<textarea id="stdFastaText" class="std-seq-text mono" spellcheck="false" placeholder=">NP1\nATGAGT...\n>NP2\n..."></textarea>
      <div class="ui-dialog-actions" style="padding:10px 0 0;"><label class="note" style="margin:0 auto 0 0;">Segment <select id="stdFastaSeg"><option>S</option><option>L</option></select>
      &nbsp;출처 <input id="stdFastaSrc" type="text" placeholder="예: KDCA 2024" style="width:140px;"></label>
      <button type="button" class="btn btn-secondary" id="stdFaCancel">취소</button><button type="button" class="btn btn-primary" id="stdFaOk">추가</button></div>`,
  });
  m.overlay.querySelector("#stdFaCancel").addEventListener("click", m.close);
  m.overlay.querySelector("#stdFaOk").addEventListener("click", () => {
    const recs = []; let cur = null;
    for (const line of m.overlay.querySelector("#stdFastaText").value.split("\n")) {
      if (line.startsWith(">")) { cur = { name: line.slice(1).trim().split(/\s+/)[0] || "", seq: "" }; recs.push(cur); }
      else if (cur) cur.seq += line.trim();
    }
    const seg = m.overlay.querySelector("#stdFastaSeg").value, src = m.overlay.querySelector("#stdFastaSrc").value.trim();
    if (!recs.length) { uiAlert("FASTA 레코드(> 로 시작하는 줄)를 찾지 못했습니다."); return; }
    const spare = stdMgr.rows.filter((r) => r.name || r.sequence);
    stdMgr.rows = spare.concat(recs.map((x) => ({ ...stdBlankRow(), name: x.name, segment: seg, sequence: stdClean(x.seq), source: src })));
    stdMgr.rows.push(stdBlankRow()); stdMgr.errors = {};
    m.close(); renderStdManager();
  });
}

async function stdSave() {
  const now = new Set(stdMgr.rows.map((r) => r.name.trim()).filter(Boolean));
  const gone = [...stdMgr.saved.keys()].filter((n) => !now.has(n));
  if (gone.length && !(await uiConfirm(`다음 표준물질의 이름이 사라집니다 (이름 변경 또는 삭제):\n\n${gone.join("\n")}\n\n이름은 결과와 실험 검증 기록의 ID라서 기존 기록과 연결이 끊어집니다. 계속 저장할까요?`, { okLabel: "저장", title: "표준물질 삭제·이름 변경" }))) {
    $("#stdSaveStatus").textContent = "저장 취소됨"; return;
  }
  $("#stdSaveStatus").textContent = "저장 중...";
  let res, data;
  try {
    res = await fetch("/api/standards/manifest", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ rows: stdMgr.rows }) });
    data = await res.json();
  } catch (e) { $("#stdSaveStatus").textContent = `저장 실패: ${e.message}`; return; }
  if (!res.ok) {
    stdMgr.errors = {}; stdMgr.general = [];
    for (const e of data.errors || []) {
      if (e.row < 0) stdMgr.general.push(e.message);
      else { (stdMgr.errors[e.row] = stdMgr.errors[e.row] || {})[e.field] = e.message; }
    }
    renderStdManager();
    $("#stdSaveStatus").textContent = `저장 실패: 입력값을 확인하세요 (${(data.errors || []).length}건)`;
    return;
  }
  stdMgr.rows = data.rows; stdMgr.rows.push(stdBlankRow());
  stdMgr.saved = new Map(data.rows.map((r) => [r.name, r.name]));
  stdMgr.errors = {}; stdMgr.general = [];
  renderStdManager();
  $("#stdSaveStatus").textContent = data.files_changed && data.files_changed.length ? "저장됨. 반영하려면 \"분석 다시하기\"를 눌러주세요." : "저장됨 (분석에 쓰이는 서열은 바뀌지 않았습니다).";
}

$("#stdMgrBtn").addEventListener("click", showStdManagerPage);

// ---------------------------------------------------------------------------------------
// UI 사용 안내 (guided tour): spotlights one control at a time and says what it is for. Started from
// the sidebar's "🧭 UI 사용 안내" button, and offered once on a first visit (localStorage). A step
// may switch the page first (show) and falls back to a centred card if its target isn't on screen.
const TOUR_KEY = "lasv_tour_done_v1";
const TOUR_STEPS = [
  { title: "LASV 프라이머/프로브 검증 대시보드", text: "라사 바이러스(LASV) 전체 게놈 데이터베이스를 기준으로 TaqMan 프라이머·프로브 세트가 얼마나 많은 게놈을 잡아내는지 보여주는 화면입니다. 1분 정도면 주요 기능을 둘러볼 수 있습니다. (←/→ 키로 이동, Esc로 종료)" },
  { target: "#overviewBtn", title: "전체 대시보드", text: "모든 프라이머 세트와 그룹의 커버리지를 한눈에 비교합니다. 리니지(LII, LIV 등)별로 얼마나 잡히는지 히트맵으로 볼 수 있고, 행을 클릭하면 그 세트의 상세로 이동합니다.", show: () => $("#overviewBtn").click() },
  { target: "#setList", title: "프라이머 세트 목록", text: "분석된 세트(이름 옆 막대 = 검증 통과율)입니다. 클릭하면 상세 화면이 열리고, ✎ 버튼으로 이름을 바꾸며, 드래그로 순서를 바꿉니다. 위 검색창으로 이름·제조사를 찾을 수 있습니다." },
  { target: "#detail .card-grid", title: "세트 상세: 요약 카드", text: "세트를 열면 Coverage, probe dropout, 표준물질 매치/대표성, 참조유전체 결합 위치 등이 카드로 나옵니다. 카드 오른쪽 위 − / ⛶ 버튼으로 접거나 전체화면으로 볼 수 있습니다.", show: async () => { const s = state.sets[0]; if (s) { await selectSet(s.set_id); await new Promise((r) => setTimeout(r, 900)); } } },
  { target: "#detail .msa-card", title: "MSA와 게놈별 결합 상세", text: "프라이머가 결합하는 구간의 다중서열정렬입니다. 게놈 행이나 아래 '게놈별 미스매치' 표의 행을 클릭하면 그 게놈에서 F/Probe/R이 어떻게 결합하는지(미스매치 위치, Tm, ΔG)를 창으로 보여줍니다." },
  { target: "#detail .tree-card", title: "Nextstrain 계통수", text: "'계통수 불러오기'를 누르면 Nextstrain 트리에 이 세트의 검증 결과를 색으로 입힙니다. 초록=통과, 빨강=불일치, 회색=우리 DB에 없음이고, 범례를 클릭하면 해당 상태만 남깁니다." },
  { target: "#addSetBtn", title: "새 프라이머 세트 분석", text: "서열(F/R/Probe)만 넣으면 전체 게놈 대비 mismatch, 표준물질 매치, 결합 위치까지 바로 계산합니다. '저장'을 체크하면 다음 전체 분석에도 포함됩니다.", show: () => $("#overviewBtn").click() },
  { target: "#primersBtn", title: "프라이머 입력 / 관리", text: "엑셀처럼 편집하는 입력표입니다. Enter/Tab으로 이동하고 엑셀에서 복사한 영역을 그대로 붙여넣을 수 있습니다. ID(P001…)는 자동 부여되어 이름을 바꿔도 유지되고, '분석진행'을 해제한 행은 분석에서 제외됩니다. 비고에는 제조사 제품명 등을 적으세요." },
  { target: "#stdMgrBtn", title: "표준물질 입력 / 관리", text: "양성대조물질(표준물질)의 서열을 등록·수정합니다. 이름이 곧 결과의 ID이고, 분석진행을 해제하면 제외됩니다. FASTA를 붙여넣어 한 번에 추가할 수도 있습니다. 저장 후 '분석 다시하기'로 반영하세요." },
  { target: "#standardsBtn", title: "표준물질 대표성", text: "KDCA 표준물질이 전체 LASV 게놈을 얼마나 대표하는지 평가합니다. 90% 이상 적합, 50% 이상 재검토, 미만이면 재선정을 권고하며, 리니지별로 어디가 비는지도 보여줍니다." },
  { target: "#diversityBtn", title: "리니지별 유전체 다양성", text: "리니지별 nucleotide diversity(π) 곡선입니다. 드래그로 확대하고 마우스를 올리면 위치별 값이 보이며, 평활화 정도를 바꿀 수 있습니다. 프라이머가 π가 낮은(보존된) 구간에 있는지 확인하는 용도입니다." },
  { target: "#metaTableBtn", title: "메타데이터 표", text: "모든 게놈의 메타데이터를 스프레드시트로 봅니다. 열 제목 클릭=정렬, ▾=값 골라 거르기(퀵 필터), 아래 입력칸=직접 필터입니다. '데이터 QC' 열과 '데이터 QC 통과만' 체크로 분석에 쓰인 게놈을 구분합니다." },
  { target: "#metadataBtn", title: "메타데이터 / 지도", text: "국가별 게놈 분포 지도와 연도별 추이, 데이터 출처·필터링 요약(몇 개 중 몇 개가 QC를 통과했는지)을 봅니다." },
  { target: "#paramsBtn", title: "분석 파라미터", text: "허용 mismatch 개수, probe 기준, 3' 말단 조건, CPU 코어 수 등을 바꿉니다. 저장 후 '분석 다시하기'를 눌러야 반영되는 항목이 있으니 화면의 안내를 확인하세요." },
  { target: "#reportBtn", title: "종합 리포트 PDF", text: "모든 세트의 결과를 세트당 1페이지로 정리한 PDF를 만듭니다." },
  { target: ".pipeline-btn-row", title: "데이터 업데이트 / 분석 다시하기", text: "'데이터 업데이트'는 Nextstrain의 최신 서열을 다시 받고, '분석 다시하기'는 현재 데이터로 전체를 다시 계산합니다. 실행 중에는 진행 로그 창이 열리며, 창을 닫아도 서버에서 계속 실행됩니다. (접속할 때 새 데이터가 있으면 자동으로 물어봅니다.)" },
  { target: "#groupModeBtn", title: "그룹 만들기", text: "같은 segment의 여러 세트를 묶어, 합집합 커버리지와 서로 다른 세트의 F/R이 교차 결합할 때의 문제까지 평가합니다." },
  { target: "#themeToggleBtn", title: "화면 설정", text: "☾/☀ 라이트·다크 모드, EN 한국어↔영어 전환이 여기 있고, 사이드바 맨 아래 '스타일' 버튼으로 Lumo 스타일과 기본 스타일을 바꿀 수 있습니다." },
  { target: "#requirementsBtn", title: "요구사항 기록", text: "프로젝트 시작부터 요청하신 사항과 검토·수정 내용, 처리 상태를 시간순으로 모아 둔 기록입니다. 검색과 구분·상태 필터로 찾아볼 수 있습니다." },
  { target: "#tourBtn", title: "다시 보려면", text: "이 안내는 언제든 이 버튼으로 다시 볼 수 있습니다. 파이프라인 자체의 설명은 바로 위 '파이프라인 개요 / 가이드'에 있습니다. 즐거운 분석 되세요!" },
];

let tourState = null;

function tourEnd(markDone = true) {
  if (!tourState) return;
  document.removeEventListener("keydown", tourState.onKey, true);
  window.removeEventListener("resize", tourState.onResize);
  tourState.root.remove();
  tourState = null;
  if (markDone) { try { localStorage.setItem(TOUR_KEY, "1"); } catch { /* private mode */ } }
}

async function tourShow(i) {
  const t = tourState;
  if (!t) return;
  t.i = Math.max(0, Math.min(TOUR_STEPS.length - 1, i));
  const step = TOUR_STEPS[t.i];
  if (step.show) { try { await step.show(); } catch { /* fall back to a centred card */ } }
  if (tourState !== t) return;
  let target = step.target ? document.querySelector(step.target) : null;
  if (target && (target.hidden || !target.getClientRects().length)) target = null;
  if (target) target.scrollIntoView({ block: "center", behavior: "instant" });
  const r = target ? target.getBoundingClientRect() : null;
  const spot = t.root.querySelector(".tour-spot");
  const card = t.root.querySelector(".tour-card");
  if (r) {
    Object.assign(spot.style, { display: "block", left: `${r.left - 6}px`, top: `${r.top - 6}px`, width: `${r.width + 12}px`, height: `${r.height + 12}px` });
  } else {
    spot.style.display = "none";
  }
  const last = t.i === TOUR_STEPS.length - 1;
  card.innerHTML = `<div class="tour-step">${t.i + 1} / ${TOUR_STEPS.length}</div>
    <div class="tour-title">${escHtml(step.title)}</div>
    <div class="tour-text">${escHtml(step.text)}</div>
    ${r ? "" : (step.target ? '<div class="note" style="margin-top:6px;">(이 항목은 지금 화면에 보이지 않습니다)</div>' : "")}
    <div class="tour-actions">
      <button type="button" class="btn btn-secondary" data-tour="skip">${last ? "닫기" : "건너뛰기"}</button>
      <span style="flex:1"></span>
      <button type="button" class="btn btn-secondary" data-tour="prev" ${t.i === 0 ? "disabled" : ""}>이전</button>
      <button type="button" class="btn btn-primary" data-tour="next">${last ? "완료" : "다음"}</button>
    </div>`;
  // place the card beside the spotlight (right of it if there is room, else below/above); centred otherwise
  const cw = Math.min(380, window.innerWidth - 24);
  card.style.width = `${cw}px`;
  let left = (window.innerWidth - cw) / 2, top = Math.max(24, window.innerHeight * 0.3);
  if (r) {
    const ch = card.offsetHeight || 200;
    if (r.right + 16 + cw <= window.innerWidth) { left = r.right + 16; top = Math.min(Math.max(12, r.top), window.innerHeight - ch - 12); }
    else if (r.bottom + 16 + ch <= window.innerHeight) { left = Math.min(Math.max(12, r.left), window.innerWidth - cw - 12); top = r.bottom + 16; }
    else { left = Math.min(Math.max(12, r.left), window.innerWidth - cw - 12); top = Math.max(12, r.top - ch - 16); }
  }
  card.style.left = `${left}px`;
  card.style.top = `${top}px`;
  card.querySelector('[data-tour="next"]').focus();
}

function startTour() {
  if (tourState) tourEnd(false);
  const root = el("div", "tour-root");
  root.innerHTML = '<div class="tour-block"></div><div class="tour-spot"></div><div class="tour-card" role="dialog" aria-modal="true"></div>';
  document.body.appendChild(root);
  const onKey = (ev) => {
    if (ev.key === "Escape") { ev.preventDefault(); ev.stopPropagation(); tourEnd(); }
    else if (ev.key === "ArrowRight") { ev.preventDefault(); ev.stopPropagation(); tourNext(); }
    else if (ev.key === "ArrowLeft") { ev.preventDefault(); ev.stopPropagation(); tourShow(tourState.i - 1); }
  };
  const onResize = () => { if (tourState) tourShow(tourState.i); };
  tourState = { root, i: 0, onKey, onResize };
  document.addEventListener("keydown", onKey, true);
  window.addEventListener("resize", onResize);
  root.addEventListener("click", (ev) => {
    const b = ev.target.closest("[data-tour]");
    if (!b || !tourState) return;
    if (b.dataset.tour === "skip") tourEnd();
    else if (b.dataset.tour === "prev") tourShow(tourState.i - 1);
    else tourNext();
  });
  tourShow(0);
}

function tourNext() {
  if (!tourState) return;
  if (tourState.i >= TOUR_STEPS.length - 1) tourEnd();
  else tourShow(tourState.i + 1);
}

$("#tourBtn").addEventListener("click", startTour);

// First visit: offer the tour once (never forced; "나중에" also counts as answered).
async function offerTourOnFirstVisit() {
  try { if (localStorage.getItem(TOUR_KEY)) return; } catch { return; }
  const go = await uiConfirm("처음 오셨나요? 화면 구성과 주요 기능을 1분 안에 안내해 드릴까요?\n(나중에는 사이드바의 '🧭 UI 사용 안내' 버튼으로 다시 볼 수 있습니다.)",
    { title: "환영합니다", okLabel: "안내 보기", cancelLabel: "나중에", danger: false });
  try { localStorage.setItem(TOUR_KEY, "1"); } catch { /* private mode */ }
  if (go) startTour();
}
document.addEventListener("click", (ev) => {
  const a = ev.target.closest && ev.target.closest(".std-rep-link");
  if (a) { ev.preventDefault(); showStandardsPage(); }
});

// One-line verdict on a set's own detail page (setData.standard_representativeness).
function renderStdRepCard(setData) {
  const r = setData.standard_representativeness;
  if (!r) return "";
  const e = r.set || {};
  let body;
  if (e.applicable) {
    const nm = (e.covering_standards || []).join(", ");
    body = `<div class="big">${stdGradeBadge(e.grade, e.pct)}</div>
      <div class="note" style="margin-top:6px;">${escHtml(nm)} 기준 결합 부위 일치 ${e.union.matched} / ${e.union.assessable}개 게놈.
      ${e.grade === "good" ? "" : e.grade === "review" ? "표준물질 재검토를 권고합니다." : "표준물질 재선정을 진지하게 고려하세요."}</div>`;
  } else {
    const u = r.segment_union || {};
    body = `<div class="note">${STD_NA_REASON[e.reason] || STD_NA_REASON.unknown}.</div>
      <div class="note" style="margin-top:6px;">참고 — 이 segment의 표준물질 전체 대표성: ${stdGradeBadge(u.grade, u.pct)}</div>`;
  }
  return `<div class="card"><h3>표준물질 대표성</h3>${body}
    <div style="margin-top:6px;"><a href="#" class="std-rep-link">자세히 보기 →</a></div></div>`;
}

// Inline footnotes for bioinformatics jargon on the guide page -- a hover-title superscript
// next to a term's FIRST mention, plus a persistent glossary list at the bottom of the page,
// so the workflow text can use precise domain terms without assuming the reader already
// knows them.
const GUIDE_FOOTNOTES = [
  { key: "iupac", term: "IUPAC 축퇴염기(degenerate base)",
    def: "한 위치에 둘 이상의 염기 중 하나가 올 수 있음을 나타내는 표준 코드(예: R=A/G, Y=C/T, N=A/C/G/T) — 리니지마다 다른 변이 위치를 프라이머 서열 하나로 함께 커버할 때 쓴다." },
  { key: "amplicon", term: "amplicon(증폭산물)",
    def: "PCR로 증폭되는 구간 전체 — Forward 프라이머 결합 위치부터 Reverse 프라이머 결합 위치까지의 서열." },
  { key: "mismatch", term: "mismatch",
    def: "프라이머/프로브 서열과 대상 게놈 서열이 서로 다른 염기의 개수. 0이면 완전일치(perfect match)." },
  { key: "tm", term: "Tm(melting temperature)",
    def: "프라이머/프로브-주형 이중나선의 절반이 풀리는 온도(융해 온도). 세트 안에서 Tm 차이가 크면 PCR 효율이 떨어질 수 있다." },
  { key: "dg", term: "&Delta;G(Gibbs 자유에너지 변화)",
    def: "값이 더 음수(더 낮음)일수록 self-dimer/cross-dimer/hairpin 같은 원치 않는 이차구조가 강하게 형성된다는 뜻." },
  { key: "pi", term: "nucleotide diversity(&#960;)",
    def: "population에서 임의의 두 서열을 뽑았을 때 그 위치의 염기가 서로 다를 평균 확률. 높을수록 해당 위치/리니지의 유전적 다양성이 크다." },
  { key: "vcf", term: "VCF(Variant Call Format)",
    def: "각 게놈이 참조유전체의 특정 위치에서 어떤 염기를 갖는지 사이트 단위로 기록하는 표준 유전체 파일 형식." },
  { key: "lineage", term: "lineage(계통, 리니지)",
    def: "계통수 상에서 공통 조상을 공유하는 게놈들의 그룹 — LASV는 LI~LVII 등 7개 주요 리니지로 분류된다." },
];
function fn(key) {
  const i = GUIDE_FOOTNOTES.findIndex((f) => f.key === key);
  if (i < 0) return "";
  return `<sup class="footnote-ref" title="${GUIDE_FOOTNOTES[i].def.replace(/&(?!amp;)/g, "&amp;").replace(/"/g, "&quot;")}">${i + 1}</sup>`;
}
function renderGuideFootnoteList() {
  return `<ol class="footnote-list">${GUIDE_FOOTNOTES.map((f) => `<li><b>${f.term}</b> — ${f.def}</li>`).join("")}</ol>`;
}

function guideStage(num, title, desc, items) {
  return `
    <div class="guide-stage">
      <div class="guide-stage-num">${num}</div>
      <div class="guide-stage-body">
        <div class="guide-stage-title">${title}</div>
        <div class="guide-stage-desc">${desc}</div>
        ${items ? `<ul class="guide-stage-items">${items.map((i) => `<li>${i}</li>`).join("")}</ul>` : ""}
      </div>
    </div>`;
}

// Process/call diagram for the guide page: WHICH script calls WHICH external tool, and how
// the two main branches (core verdict path vs. dashboard-only supplementary analyses) relate
// -- the textual step list above already explains WHAT each stage computes, this shows HOW
// the actual process calls wire together (Snakemake DAG dependency + external binary calls).
const DIAG_TOOL_COLOR = { mfeprimer: "#F2A65A", mafft: "#54A24B", pixy: "#B279A2", python: "var(--muted)" };
// The diagram texts are written with HTML entities (&rarr; ...), so only a bare "&" and "<" / ">" are escaped.
function diagText(t) {
  return String(t).replace(/&(?![a-zA-Z]+;|#\d+;)/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
}
function diagBox(x, y, w, h, lines, opts) {
  opts = opts || {};
  const fill = opts.fill || "var(--card)";
  const stroke = opts.stroke || "var(--border)";
  const dashed = opts.dashed ? ' stroke-dasharray="5,3"' : "";
  const textColor = opts.textColor || "var(--fg)";
  const lineH = opts.lineH || 13;
  const fontSize = opts.fontSize || 11;
  const startY = y + h / 2 - ((lines.length - 1) * lineH) / 2 + fontSize / 3;
  const textLines = lines.map((line, i) =>
    `<text x="${x + w / 2}" y="${(startY + i * lineH).toFixed(1)}" text-anchor="middle" font-size="${line.small ? fontSize - 1.5 : fontSize}" font-weight="${line.bold ? 700 : 400}" fill="${line.color || textColor}">${diagText(line.text || line)}</text>`
  ).join("");
  return `<g>
    <rect x="${x}" y="${y}" width="${w}" height="${h}" rx="7" fill="${fill}" stroke="${stroke}" stroke-width="1.2"${dashed}></rect>
    ${textLines}
  </g>`;
}
function diagArrow(x1, y1, x2, y2, opts) {
  opts = opts || {};
  const stroke = opts.stroke || "var(--muted)";
  const dashed = opts.dashed ? ' stroke-dasharray="4,3"' : "";
  const midX = (x1 + x2) / 2, midY = (y1 + y2) / 2;
  const angle = Math.atan2(y2 - y1, x2 - x1);
  const ah = 5;
  const p1x = x2 - ah * Math.cos(angle - Math.PI / 7), p1y = y2 - ah * Math.sin(angle - Math.PI / 7);
  const p2x = x2 - ah * Math.cos(angle + Math.PI / 7), p2y = y2 - ah * Math.sin(angle + Math.PI / 7);
  return `<line x1="${x1}" y1="${y1}" x2="${x2}" y2="${y2}" stroke="${stroke}" stroke-width="1.4"${dashed}></line>
    <polygon points="${x2},${y2} ${p1x.toFixed(1)},${p1y.toFixed(1)} ${p2x.toFixed(1)},${p2y.toFixed(1)}" fill="${stroke}"></polygon>
    ${opts.label ? `<text x="${midX + (opts.labelDx || 6)}" y="${midY + (opts.labelDy || -4)}" font-size="9.5" fill="${stroke}">${escHtml(opts.label)}</text>` : ""}`;
}
function diagToolPill(x, y, tool, label) {
  const color = DIAG_TOOL_COLOR[tool] || "var(--muted)";
  return `<text x="${x}" y="${y}" text-anchor="middle" font-size="9.5" font-weight="700" fill="${color}">&#9656; ${diagText(label)}</text>`;
}

function renderPipelineDiagram() {
  const W = 1000;
  const parts = [];
  const B = (x, y, w, h, lines, opts) => parts.push(diagBox(x, y, w, h, lines, opts));
  const A = (x1, y1, x2, y2, opts) => parts.push(diagArrow(x1, y1, x2, y2, opts));

  // ---- vertical spine: shared ingestion, common to every branch below ----
  B(390, 10, 220, 44, [{ text: "Nextstrain (nextstrain.org/lassa)", bold: true }, { text: "S/L segment 서열 + lineage 메타데이터", small: true }],
    { fill: "var(--card)", dashed: true });
  A(500, 54, 500, 74, {});
  B(390, 76, 220, 40, [{ text: "fetch_raw.sh", bold: true }, { text: "curl로 원본 다운로드", small: true }]);
  A(500, 116, 500, 136, {});
  B(360, 138, 280, 40, [{ text: "rename_headers.py &rarr; filter.py (+ 방향 통일)", bold: true }, { text: "N 염기 / 길이 기준 필터링", small: true }]);
  A(500, 178, 500, 200, {});
  B(340, 202, 320, 44, [{ text: "population FASTA (S/L, filtered) + mfeprimer index", bold: true }],
    { fill: "var(--card)", dashed: true, textColor: "var(--accent)" });
  A(500, 246, 500, 268, { label: "이 지점에서 두 경로로 분기", labelDx: 10, labelDy: 4 });

  // ---- fan-out bus line ----
  parts.push(`<line x1="150" y1="268" x2="850" y2="268" stroke="var(--border)" stroke-width="1.2"></line>`);
  const branchXs = [150, 850];
  branchXs.forEach((bx) => A(500, 268, bx, 268, { stroke: "var(--border)" }));

  // ---- LEFT column: core verdict path ----
  const lx = 40, lw = 420;
  parts.push(`<text x="${lx + lw - 4}" y="292" text-anchor="end" font-size="11" font-weight="700" fill="var(--accent)">핵심 판정 경로 (reconcile까지 이어짐)</text>`);
  A(150, 268, 150, 300, {});
  B(lx, 300, lw, 52, [{ text: "primer_match", bold: true }, { text: "run_primer.py — mfeprimer --virus (mismatch ≤2)", small: true }], { stroke: "var(--accent)" });
  parts.push(diagToolPill(lx + lw - 70, 316, "mfeprimer", "mfeprimer --virus"));
  A(150, 352, 150, 384, {});
  B(lx, 384, lw, 60, [{ text: "build_alignment_table.py", bold: true }, { text: "두 엔진(mfeprimer + 게놈 직접 탐색) 비교, 게놈 염기로 mismatch 재계산", small: true }], { stroke: "var(--accent)" });
  parts.push(diagToolPill(lx + lw - 70, 400, "mfeprimer", "mfeprimer search"));
  parts.push(diagToolPill(lx + lw - 70, 412, "python", "+ 자체 Python"));
  A(150, 444, 150, 476, {});
  B(lx, 476, lw, 44, [{ text: "reconcile.py", bold: true }, { text: "validated = primer(3′ 끝 포함) &cap; probe(amplicon 안)", small: true }], { stroke: "var(--accent)" });
  A(150, 520, 150, 552, {});
  B(lx, 552, lw, 44, [{ text: "aggregate.py / aggregate_groups.py", bold: true }, { text: "세트·그룹 전체 요약(summary.tsv 등)", small: true }], { stroke: "var(--accent)" });

  // ---- RIGHT column: supplementary, dashboard-only analyses ----
  const rx = 540, rw = 420;
  parts.push(`<text x="${rx + 4}" y="292" text-anchor="start" font-size="11" font-weight="700" fill="var(--muted)">부가 분석 (대시보드 전용 &mdash; reconcile을 거치지 않음)</text>`);
  const rRows = [
    { y: 300, lines: [{ text: "refgenome_binding.py", bold: true }, { text: "단일 참조유전체 결합위치", small: true }], tool: "mfeprimer", toolLabel: "mfeprimer --virus" },
    { y: 358, lines: [{ text: "standard_primer/probe_match &rarr; standard_alignment_table.py", bold: true }, { text: "KDCA 표준물질 DB 검증", small: true }], tool: "mfeprimer", toolLabel: "mfeprimer" },
    { y: 416, lines: [{ text: "msa_alignment.py (build_msa.py)", bold: true }, { text: "세트별 amplicon 다중정렬", small: true }], tool: "mafft", toolLabel: "mafft" },
    { y: 474, lines: [{ text: "dimer_hairpin_check.py", bold: true }, { text: "self/cross-dimer, hairpin 위험", small: true }], tool: "mfeprimer", toolLabel: "mfeprimer dimer+hairpin" },
    { y: 532, lines: [{ text: "build_diversity_vcf.py &rarr; run_pixy.py", bold: true }, { text: "전체게놈 정렬(방향보정) &rarr; 리니지별 &#960;", small: true }], tool: "pixy", toolLabel: "mafft &rarr; pixy" },
  ];
  rRows.forEach((r) => {
    A(850, 268, 850, r.y, {});
    B(rx, r.y, rw, 50, r.lines, { stroke: "var(--border)" });
    parts.push(diagToolPill(rx + rw - 78, r.y + 16, r.tool, r.toolLabel));
  });

  // ---- convergence into shared outputs ----
  const convY = 616;
  A(150, 596, 150, convY, {});
  rRows.forEach((r) => A(850, r.y + 50, 850, convY, {}));
  parts.push(`<line x1="150" y1="${convY}" x2="850" y2="${convY}" stroke="var(--border)" stroke-width="1.2"></line>`);
  A(500, convY, 500, convY + 22, {});
  B(300, convY + 24, 400, 44, [{ text: "output/summary/*, output/dashboard/*", bold: true }, { text: "모든 단계의 최종 산출물", small: true }],
    { fill: "var(--card)", dashed: true, textColor: "var(--accent)" });
  A(500, convY + 68, 500, convY + 90, {});
  B(260, convY + 92, 480, 50, [{ text: "dashboard.py (정적 SVG) &nbsp;&Vert;&nbsp; webapp/server.py (인터랙티브)", bold: true }], { stroke: "var(--accent)", fill: "rgba(76,120,168,.10)" });

  // ---- side note: on-demand path that bypasses Snakemake entirely ----
  const noteY = convY + 160;
  B(60, noteY, 880, 56,
    [{ text: "온디맨드 예외: 웹앱의 \"+ 새 프라이머 세트 분석\", 그룹 이차구조/미세조정 제안", bold: true },
     { text: "Snakemake DAG를 거치지 않고 live_analyze.py / webapp/server.py가 그 자리에서 직접 mfeprimer(--virus/search/dimer/hairpin) + mafft를 호출", small: true }],
    { stroke: "var(--warn)", dashed: true });

  const H = noteY + 76;
  return `<svg viewBox="0 0 ${W} ${H}" class="pipeline-diagram-svg" preserveAspectRatio="xMidYMid meet">${parts.join("")}</svg>
    <div class="pipeline-diagram-legend">
      <span><span class="dot" style="background:${DIAG_TOOL_COLOR.mfeprimer}"></span>mfeprimer 호출</span>
      <span><span class="dot" style="background:${DIAG_TOOL_COLOR.mafft}"></span>mafft 호출</span>
      <span><span class="dot" style="background:${DIAG_TOOL_COLOR.pixy}"></span>pixy 호출</span>
      <span><span class="dot" style="background:var(--muted)"></span>순수 Python(외부 도구 미호출)</span>
    </div>`;
}

function renderGuidePage() {
  const flow = [
    guideStage("1", "데이터 수집 (fetch)",
      "Nextstrain(nextstrain.org/lassa)에서 LASV S segment / L segment의 최신 서열과 메타데이터(lineage, 국가, 날짜 등)를 내려받습니다."),
    guideStage("2", "필터링 (filter)",
      "모호 염기(N)가 하나라도 있거나 길이가 기준(S&ge;3000bp, L&ge;6000bp) 미만인 게놈은 제외해, mfeprimer 검색에 쓰이는 population DB를 만듭니다. 통과한 게놈은 <b>참조유전체와 같은 방향으로 통일</b>합니다 — 반대 가닥으로 등록된 서열은 참조와 공유하는 11-mer 개수로 판별해 역상보로 뒤집은 뒤 분석하므로, 좌표·증폭 산물 서열·정렬이 모든 게놈에서 같은 방향입니다(뒤집은 게놈 목록은 <code>data/filtered/&lt;seg&gt;_orientation.tsv</code>)."),
    guideStage("3~4", "F/R · 프로브 결합 탐색 (primer_match / probe_match / alignment_table)",
      `각 세트의 Forward/Reverse 프라이머 쌍은 <b>두 가지 엔진</b>으로 모든 게놈에서 찾습니다: <b>mfeprimer --virus</b>(mismatch&le;2로 넓게 탐색)와 <b>게놈 전체 직접 탐색</b>(모든 위치·양쪽 가닥을 직접 비교). Probe도 같은 방식으로 mfeprimer <b>search</b>와 직접 탐색을 함께 씁니다. 두 결과를 비교해 더 나은 쪽을 쓰고, 어느 쪽 결과인지 기록합니다.`,
      [`mismatch${fn("mismatch")}는 mfeprimer가 보고한 서열이 아니라 <b>그 위치의 실제 게놈 염기</b>와 프라이머를 비교해 셉니다. IUPAC 축퇴염기${fn("iupac")}(R/Y/S/W/K/M/B/D/H/V/N)는 어느 쪽 염기든 해당 코드에 포함되면 일치로 봅니다.`,
       `F/R 프라이머 쌍은 두 프라이머가 서로 마주 보는 유효한 amplicon${fn("amplicon")}(방향·최대 2000 bp)을 이루는 조합만 인정하고, 프라이머마다 <b>3′ 끝 연속 일치</b>(기본 2염기)도 확인합니다. Probe는 <b>두 프라이머 사이(amplicon 안)</b>에 있어야 하며 양쪽 가닥을 모두 봅니다.`,
       `mismatch 0/1/2 단계별 커버리지는 따로 검색하지 않고, 위 게놈별 표(<code>mismatch_table.tsv</code>)에서 계산합니다.`,
       `비교 결과는 게놈별 <code>fr_source</code> 열, <code>engine_comparison.tsv</code>, 세트별 일치율 <code>output/summary/engine_concordance.tsv</code>에 남습니다.`]),
    guideStage("5", "판정 (reconcile / aggregate)",
      "<b>validated = (F/R 프라이머 mismatch &le; 설정값, 3′ 끝 일치) AND (probe mismatch = 설정값, amplicon 안)</b> 조건을 만족하는 게놈만 최종 통과로 집계합니다. 설정값은 \"⚙ 분석 파라미터 설정\" 페이지에서 바꿀 수 있고, 바꾼 뒤 \"분석 다시하기\"를 누르면 mfeprimer를 다시 돌리지 않고 판정 단계만 다시 계산합니다.",
      ["프라이머만 통과하고 프로브가 실패한 게놈은 <b>probe dropout</b>으로 따로 셉니다(프로브가 맞지 않으면 프라이머가 증폭해도 신호가 안 날 수 있기 때문).",
       "세트별 요약(<code>summary.tsv</code>)과 리니지별 요약(<code>summary_by_lineage.tsv</code>)이 만들어집니다."]),
    guideStage("5b", "판정 방식과 mfeprimer의 관계 (검증 결과)",
      "mfeprimer는 <b>후보 결합 위치를 찾는 데</b> 쓰고, 판정에 쓰는 mismatch 개수는 그 위치의 <b>실제 게놈 염기와 프라이머를 직접 비교해 다시 셉니다</b>. mfeprimer가 보고하는 mismatch만으로는 5′ 쪽 불일치를 놓치기 때문입니다.",
      ["<b>왜 다시 세는가</b>: mfeprimer의 <code>--misMatch</code>는 프라이머 전체가 아니라 <b>3′ 끝 k-mer(기본 9염기) 안</b>의 mismatch만 반영하고, 보고하는 서열(fpSeq/rpSeq)도 게놈이 아니라 프라이머 자신입니다. 그래서 5′ 쪽에 mismatch가 있어도 \"0\"으로 나옵니다. 예) Nipah P001의 R 프라이머는 78개 게놈에서 5′ 끝 두 번째 염기가 다릅니다(mfeprimer는 0으로 보고).",
       "<b>검증 방법</b>: 게놈의 모든 위치·양쪽 가닥에 프라이머를 한 칸씩 대어 mismatch를 세는 독립 계산(<code>scripts/audit_recount.py</code>)과 파이프라인의 게놈별 표를 비교합니다. LASV 25개 세트와 Nipah 2개 세트의 모든 게놈에서 표가 독립 계산보다 적은 mismatch를 주장하는 경우가 없습니다(표가 더 높게 나오는 경우는 두 프라이머가 증폭 산물을 이루는 사이트만 고르기 때문이며 정상).",
       "<b>mfeprimer 원 결과와의 비교(2026-10-15)</b>: ① mfeprimer만 통과시키는 게놈은 5′ 쪽 mismatch가 2개 이상인 경우입니다(LASV P003 98개, P011 208개). ② 우리만 통과시키는 게놈은 mfeprimer가 <code>--misMatch 1</code>에서 실제로 mismatch 1개 이하인 사이트를 놓친 경우입니다(LASV P001 90개). ③ <code>--misMatch 2</code>는 거의 모든 게놈을 통과시켜 판정 기준으로 쓸 수 없습니다(LASV P016~P019는 우리 기준 0건인데 mfeprimer는 100%). ④ 프로브 단독 검색은 두 방법이 거의 같고, 차이(LASV P007 88→39)는 증폭 산물 밖에서 결합하는 프로브를 우리가 제외한 결과입니다.",
       "<b>두 엔진을 모두 돌려 비교</b>: 모든 게놈에 mfeprimer와 <b>게놈 전체 직접 탐색</b>(모든 위치·양쪽 가닥, 마주 보는 쌍, 2000 bp 이내)을 함께 적용하고, 더 나은 쌍(3′ 끝 정상 → 총 mismatch 적은 쪽)을 최종 결과로 씁니다. 어느 엔진의 쌍인지는 게놈별 표의 <code>fr_source</code>에, 게놈마다 두 엔진의 결과와 판정 일치 여부는 <code>output/primers/&lt;ID&gt;/alignment/engine_comparison.tsv</code>와 요약 <code>output/summary/engine_concordance.tsv</code>에 남습니다. 직접 탐색이 mfeprimer보다 나쁜 경우는 한 건도 없었고(LASV 25개·Nipah 2개 세트), mfeprimer가 놓친 쌍은 일부 세트에서만 나옵니다(예: LASV P004 29개, P020 22개 게놈).",
       "<b>프로브 위치</b>: 프로브는 두 프라이머 사이(증폭 산물 안)에 결합해야 통과로 셉니다.",
       "<b>한계</b>: Tm·&Delta;G·dimer·hairpin은 mfeprimer의 계산값을 그대로 보여 줍니다. 이 값은 3′ 끝 변이형 기준이라 5′ 쪽 mismatch가 있는 게놈에서는 근사값입니다. 판정(validated)에는 쓰이지 않습니다.",
       "인용: Wang K, et al. MFEprimer-3.0: quality control for PCR primers. <i>Nucleic Acids Res</i>. 2019 (doi:10.1093/nar/gkz351)."]),
    guideStage("6", "부가 분석 (같은 탐색 로직 재사용)",
      "판정과 별개로, 검증에 도움이 되는 추가 정보를 같은 방식으로 계산합니다.",
      ["<b>표준물질 매치(선택 사항)</b>: 분석 대상 게놈 대신 표준물질(양성대조물질, 몇 개의 짧은 서열) DB에 대해 동일하게 검색 — 실제 wet-lab 양성대조군과 맞는지 확인. 표준물질을 등록하지 않으면 이 단계는 건너뛰고 화면에는 \"미등록\"으로 표시됩니다.",
       `<b>참조유전체 결합위치</b>: 단일 대표 유전체(${typeof refAccList === "function" ? refAccList() : "S=AY628203, L=AY628202"})에 대해 F/Probe/R이 어느 유전자 영역에 결합하는지 좌표로 확인.`,
       `<b>게놈별 F/Probe/R 미스매치 표</b>: 분석 대상 게놈 각각에 대해 F/Probe/R mismatch, 결합위치, Tm${fn("tm")}, &Delta;G${fn("dg")}(결합 자유에너지)를 계산(Tm·&Delta;G는 mfeprimer 값).`,
       "<b>전체 게놈 정렬(MSA)</b>: mafft로 모든 결합 게놈의 amplicon을 한 번에 정렬. 점(.)은 프라이머/프로브 구간에서는 올리고 자신의 염기, 그 밖에서는 consensus와의 일치를 뜻하며 선택 상자로 consensus 기준으로 바꿀 수 있습니다.",
       "<b>그룹 교차 결합 분석</b>: 여러 세트를 그룹으로 묶었을 때, 각 세트를 독립적으로 보는 대신 그룹의 모든 F/R/Probe 올리고를 서로 짝지어(예: 세트1의 F + 세트2의 R + 어느 세트의 프로브) 두 엔진으로 모든 게놈에서 찾습니다. 그룹 커버리지는 <b>합집합(세트별 단독)</b>과 <b>교차 결합 포함</b> 두 값으로 함께 표시합니다.",
       "<b>유전자 보존영역</b>: 아래 7단계의 자리별 &#960;로 참조유전체의 유전자마다 전체와 리니지/클레이드별 보존 정도, 가장 긴 보존 구간을 계산(새 분석 단계 없음)."]),
    guideStage("7", `리니지${fn("lineage")}별 유전체 다양성 (pixy &#960;)`,
      `개별 세트의 mfeprimer 매치와는 별개로, population 게놈 하나하나를 <b>mafft --keeplength --add</b>로 참조유전체 하나에 독립적으로 정렬한 뒤(LASV S/L 게놈은 GenBank 등록 방향이 서로 달라 --adjustdirection도 함께 사용) 참조 행을 기준으로 all-sites VCF${fn("vcf")}를 직접 만들어 <a href="https://github.com/ksamuk/pixy" target="_blank" rel="noopener">pixy</a>로 리니지별 nucleotide diversity(&#960;)${fn("pi")}를 1bp(염기 단위) 슬라이딩 윈도우로 계산합니다. 전체 population을 한 번에 조인트 정렬하지 않는 이유: 리니지 간 발산도가 큰 게놈이 섞이면 정렬 자체가 뒤틀려 실제로는 정상인 게놈도 결측(N) 처리되는 문제가 있었음(실측 확인 후 이 방식으로 교체).`,
      ["S/L segment 각각 하나의 all-sites VCF + populations.txt를 만들고(<code>build_diversity_vcf</code>), pixy가 리니지별 &#960;를 윈도우 단위로 산출합니다(<code>run_pixy</code>, multiallelic 사이트도 <code>--include_multiallelic_snps</code>로 포함).",
       "결과는 \"🧬 리니지별 유전체 다양성\" 페이지에서 Nextstrain의 genome-wide entropy 패널과 비슷한 형태(리니지별 색상 라인 + 유전자 경계)로 볼 수 있고, 각 프라이머 세트 상세 페이지에서도 그 세트의 F/Probe/R 결합 위치에 맞춰 확대된 오버레이 차트와 전체 게놈 정렬(MSA) 위 막대그래프(히스토그램)로 함께 표시됩니다.",
       "윈도우를 1bp(염기 단위)로 잡은 이유: amplicon 크기가 보통 100~170bp 정도라, 50bp 단위로는 한 amplicon 안에 막대가 2~3개밖에 안 나와 구조를 보기 어렵기 때문입니다. 대신 염기 단위라 리니지별 오버레이 차트는 실제 변이가 있는 정확한 위치에서만 뾰족하게 튀는 형태가 됩니다(구간 신호가 아니라 그 염기 하나의 변이 여부 신호이기 때문) — MSA 위 히스토그램은 리니지 평균을 내기 때문에 상대적으로 덜 들쭉날쭉합니다."]),
    guideStage("8", "대시보드 출력과 입력/관리",
      "위 모든 결과를 이 인터랙티브 webapp(현재 화면)에서 세트별/그룹별로 드릴다운하며 볼 수 있게 합니다. 정적 SVG 대시보드(output/dashboard/index.html)도 별도로 생성됩니다.",
      ["세트 상세 페이지: 서열/커버리지/표준물질 매치/참조유전체 결합위치/리니지 분포/전체 게놈 정렬(MSA, &#960; 히스토그램 포함)/계통수 대비 커버리지를 한 화면에서 확인. 표준물질 카드와 참조유전체 결합위치 카드는 눌러서 상세를 모달로 볼 수 있습니다.",
       "<b>세트별·그룹별 분석/재분석 버튼</b>: 새 프라이머를 올린 뒤 전체를 돌리지 않고 그 세트(또는 그룹)만 분석하거나, 입력이 같아도 강제로 다시 계산할 수 있습니다.",
       "<b>프라이머/표준물질 입력 페이지</b>: 엑셀처럼 편집하고, 엑셀(.xlsx) 파일로 내려받거나 불러올 수 있습니다(불러오면 표에 채우기만 하고 저장은 직접 누릅니다).",
       "형광물질(채널) 기록과 출처(citation, 예: 공개 논문 기반 프라이머의 MLA 인용) 표기를 세트별로 저장.",
       "실험(wet-lab) 검증 기록: 실제 RT-PCR 결과(검출 여부, Ct, 형광물질, 이미지)를 인실리코 예측과 비교.",
       "\"🗺 메타데이터 / 지도\" 페이지: 분석 대상 데이터의 수집/필터링 통계와, 국가별 보고 유전체 수를 지도로(국가에 마우스를 올리면 연도별 추이 그래프 표시).",
       "\"📄 종합 리포트 PDF 생성\": 위 내용을 워크플로 설명 + 전체 요약 + 세트당 1페이지로 정리한 PDF로 내보내기(서열은 자동으로 가려서 출력).",
       "점검: <code>./test.sh</code>(핵심 기능을 정답을 아는 합성 게놈으로 검증)와 <code>./test.sh --full</code>(실제 결과 감사: 프라이머 mismatch 재계산, 쌍 누락, 가닥 방향)."]),
  ];
  const flowHtml = flow.map((stage, i) => stage + (i < flow.length - 1 ? '<div class="guide-arrow">↓</div>' : "")).join("");

  return `
    <h2>파이프라인 개요 / 가이드</h2>
    <div class="subtitle">이 도구가 무엇을 어떻게 계산하는지, 파라미터와 결과를 어떻게 해석하면 되는지 정리한 페이지입니다.</div>

    <div class="card" style="margin-bottom:20px;">
      <h3>이 도구는 무엇인가요?</h3>
      <p>LASV(라사바이러스) TaqMan primer/probe 세트가 Nextstrain의 실제 population 유전체 데이터에 대해
      얼마나 잘 매치되는지를 자동으로 검증하는 Snakemake 파이프라인 + 인터랙티브 대시보드입니다.
      설계한 프라이머/프로브가 실제 유전적 다양성(lineage별 변이 등)을 얼마나 커버하는지, 어떤 게놈에서
      mismatch가 발생하는지, 표준물질과는 맞는지, 여러 세트를 함께 쓸 때 의도치 않게 서로 결합하지는
      않는지를 한 곳에서 확인할 수 있도록 만들었습니다.</p>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>NCBI Primer-BLAST와 무엇이 다른가요?</h3>
      <p class="note" style="margin-top:0;">둘 다 "mismatch를 검사한다"는 점은 같지만, 검증 방향이 반대입니다:
      Primer-BLAST는 <b>오프타겟(off-target) 특이성</b>(엉뚱한 곳에 안 붙는지)을 보고, 이 도구는
      <b>population 커버리지</b>(표적 바이러스의 실제 리니지/변이체에 빠짐없이 붙는지)를 봅니다.</p>
      <table class="mini-reference-table">
        <tr><th></th><th>NCBI Primer-BLAST</th><th>이 파이프라인</th></tr>
        <tr><td>목적</td><td>오프타겟 특이성 검사(엉뚱한 유전자에 안 붙는지)</td><td>population 커버리지 검증(표적 바이러스 변이체에 다 붙는지)</td></tr>
        <tr><td>비교 대상</td><td>RefSeq 등 genome-wide 데이터베이스</td><td>큐레이션된 population 게놈(S 662개 / L 366개, 7개 리니지)</td></tr>
        <tr><td>탐색 알고리즘</td><td>BLASTN(word_size=7) + Needleman-Wunsch — mismatch가 퍼져 있어 연속 일치 구간이 7bp 미만이면 놓칠 수 있음</td><td>mfeprimer --virus(바이러스 검색 특화) + 게놈 전체 직접 탐색, 두 결과를 비교</td></tr>
        <tr><td>IUPAC 축퇴염기${fn("iupac")}</td><td>축퇴 프라이머를 다루도록 설계되지 않음</td><td>IUPAC 코드를 인식해 게놈 염기와 직접 비교(축퇴 조합을 따로 확장하지 않아도 정확)</td></tr>
        <tr><td>TaqMan 프로브</td><td>별도 개념 없음(F/R 쌍 설계·검사 중심)</td><td>프로브 완전일치를 F/R과 분리된 3차 판정 기준으로 명시</td></tr>
        <tr><td>출력 단위</td><td>단일 설계/검사에 대한 후보 리스트 + 특이성 리포트</td><td>population 개별 게놈 결과 + 리니지별 % + 계통수 + 다양성(&#960;) 오버레이</td></tr>
      </table>
      <div class="note" style="margin-top:8px;">
        Primer-BLAST의 word_size=7 seed 한계 예시: 20bp 프라이머에서 mismatch가 7번·14번 위치처럼 퍼져 있으면
        (연속 일치 구간이 전부 6bp 이하가 되면) mismatch 2개뿐이어도 그 타겟을 통째로 놓칩니다 — SNIPPY/minimap2에서
        확인했던 seed 기반 aligner의 한계(12장 참고)와 같은 종류의 문제입니다.
        (<a href="https://www.idtdna.com/page/support-and-education/decoded-plus/tips-for-using-blast-to-locate-pcr-primers" target="_blank" rel="noopener">VendorD: How to Use Primer-BLAST</a>,
        <a href="https://link.springer.com/article/10.1186/1471-2105-13-134" target="_blank" rel="noopener">Ye et al. 2012, BMC Bioinformatics</a>)
      </div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>전체 워크플로우</h3>
      <div class="guide-flow">${flowHtml}</div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>프로세스 호출 다이어그램</h3>
      <p class="note" style="margin-top:0;">위 워크플로우가 "무엇을 계산하는지"였다면, 이 다이어그램은 "어느 스크립트가 실제로 어떤
        외부 프로그램을 호출하는지"를 보여줍니다. Nextstrain에서 받은 population 서열이 필터링을 거친 뒤,
        <b style="color:var(--accent);">핵심 판정 경로</b>(mfeprimer로 F/R을 찾고 → probe를 재확인하고
        → 판정하고 → 요약하는, <code>validated_pct</code>에 직접 영향을 주는 한 줄기)와
        <b>부가 분석</b>(참조유전체 결합위치, 표준물질 매치, 전체 게놈 정렬, 이차구조, 리니지별
        다양성 — 전부 대시보드 표시용이고 <code>reconcile</code>의 판정 자체에는 들어가지 않음)으로
        나뉘어 각각 mfeprimer/mafft/pixy를 호출한 뒤, 두 경로 모두 최종 결과물로 합쳐집니다. 맨 아래
        점선 박스는 Snakemake DAG를 거치지 않고 웹앱이 그 자리에서 즉시 같은 도구들을 호출하는
        예외 경로입니다("+ 새 프라이머 세트 분석", 그룹 이차구조/미세조정 제안 등).</p>
      <div class="pipeline-diagram-wrap">${renderPipelineDiagram()}</div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>사용된 프로그램</h3>
      <div class="guide-tools">
        <div class="guide-tool"><div class="guide-tool-name">Snakemake</div><div class="guide-tool-desc">전체 분석 단계를 의존성 DAG로 관리하는 워크플로우 엔진. 입력 파일이 바뀐 단계만 다시 계산합니다.</div></div>
        <div class="guide-tool"><div class="guide-tool-name">mfeprimer (v4.2.4)</div><div class="guide-tool-desc">프라이머 쌍/프로브가 주어진 서열 DB의 어디에 몇 개의 mismatch로 결합하는지 찾는 핵심 엔진. Tm(결합 온도)과 ΔG(결합 자유에너지)도 함께 계산합니다.</div></div>
        <div class="guide-tool"><div class="guide-tool-name">mafft</div><div class="guide-tool-desc">여러 게놈의 amplicon 서열을 한 번에 정렬(multiple sequence alignment)하는 도구 — "전체 게놈 정렬" 뷰와, population 전체 게놈을 방향까지 맞춰 정렬하는 리니지별 다양성 분석(--adjustdirection)에 모두 사용.</div></div>
        <div class="guide-tool"><div class="guide-tool-name">pixy</div><div class="guide-tool-desc">결측치(gap/N)가 있는 population 데이터에서도 편향 없이 nucleotide diversity(&#960;)를 추정하는 도구. 리니지별/윈도우별 &#960; 계산에 사용 — "🧬 리니지별 유전체 다양성" 페이지와 세트 상세 페이지의 &#960; 오버레이의 데이터 출처.</div></div>
        <div class="guide-tool"><div class="guide-tool-name">Nextstrain</div><div class="guide-tool-desc">LASV S/L segment 공개 서열 + 계통(lineage) 메타데이터 및 계통수의 출처.</div></div>
      </div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>레퍼런스 / 참고문헌 (도구 버전 &amp; 인용)</h3>
      <div class="note" style="margin-bottom:10px;">이 파이프라인이 실제로 사용 중인 각 도구의 버전과, 학술적으로 인용할 때 쓰는 출처입니다
        (이 환경에서 직접 확인한 버전 — <code>--version</code>/<code>version</code> 명령 결과).</div>
      <table class="mini-reference-table">
        <thead><tr><th>도구</th><th>버전</th><th>용도</th><th>인용/출처</th></tr></thead>
        <tbody>
          <tr>
            <td>mfeprimer</td><td class="mono">v4.2.4</td>
            <td>F/R/Probe 결합 위치·mismatch·Tm·ΔG 탐색, self/cross-dimer 및 hairpin 검사 (핵심 엔진)</td>
            <td>Wang, K., Li, H., Xu, Y., Shao, Q., Yi, J., Wang, R., … Qu, W. (2019). MFEprimer-3.0:
              quality control for PCR primers. <i>Nucleic Acids Research</i>.
              <a href="https://doi.org/10.1093/nar/gkz351" target="_blank" rel="noopener">doi:10.1093/nar/gkz351</a></td>
          </tr>
          <tr>
            <td>mafft</td><td class="mono">v7.526 (2024/Apr/26)</td>
            <td>Amplicon MSA("전체 게놈 정렬" 뷰), population 전체 게놈 방향 정규화(<code>--adjustdirection</code>,
              리니지별 다양성 분석용)</td>
            <td>Katoh, K., &amp; Standley, D. M. (2013). MAFFT multiple sequence alignment software
              version 7: improvements in performance and usability. <i>Molecular Biology and Evolution</i>,
              30(4), 772&ndash;780.
              <a href="https://doi.org/10.1093/molbev/mst010" target="_blank" rel="noopener">doi:10.1093/molbev/mst010</a></td>
          </tr>
          <tr>
            <td>pixy</td><td class="mono">v2.2.3</td>
            <td>리니지별 nucleotide diversity(&#960;) 추정 ("🧬 리니지별 유전체 다양성" 페이지, &#960; 오버레이,
              미세조정 제안의 데이터 출처)</td>
            <td>Korunes, K. L., &amp; Samuk, K. (2021). pixy: Unbiased estimation of nucleotide diversity
              and divergence in the presence of missing data. <i>Molecular Ecology Resources</i>, 21(4),
              1359&ndash;1368.
              <a href="https://doi.org/10.1111/1755-0998.13326" target="_blank" rel="noopener">doi:10.1111/1755-0998.13326</a></td>
          </tr>
          <tr>
            <td>htslib (bgzip/tabix)</td><td class="mono">v1.23.1</td>
            <td>pixy 입력용 VCF 압축/색인</td>
            <td>Bonfield, J. K., et al. (2021). HTSlib: C library for reading/writing high-throughput
              sequencing data. <i>GigaScience</i>, 10(2), giab007.
              <a href="https://doi.org/10.1093/gigascience/giab007" target="_blank" rel="noopener">doi:10.1093/gigascience/giab007</a></td>
          </tr>
          <tr>
            <td>Snakemake</td><td class="mono">v9.27.0</td>
            <td>전체 파이프라인 워크플로우 관리(DAG 기반 의존성/재실행)</td>
            <td>M&ouml;lder, F., et al. (2021). Sustainable data analysis with Snakemake.
              <i>F1000Research</i>, 10, 33.
              <a href="https://doi.org/10.12688/f1000research.29032.2" target="_blank" rel="noopener">doi:10.12688/f1000research.29032.2</a></td>
          </tr>
          <tr>
            <td>d3.js</td><td class="mono">v7.9.0</td>
            <td>인터랙티브 차트(참조유전체 지도, 연도별 추이, 다양성 &#960; 확대/토글 등)</td>
            <td>Bostock, M., Ogievetsky, V., &amp; Heer, J. (2011). D&sup3; Data-Driven Documents.
              <i>IEEE Transactions on Visualization and Computer Graphics</i>, 17(12), 2301&ndash;2309.
              <a href="https://doi.org/10.1109/TVCG.2011.185" target="_blank" rel="noopener">doi:10.1109/TVCG.2011.185</a> &middot;
              <a href="https://d3js.org" target="_blank" rel="noopener">d3js.org</a></td>
          </tr>
          <tr>
            <td>Python</td><td class="mono">3.12.10</td>
            <td>파이프라인 스크립트 + 웹앱 서버 전체</td>
            <td>Van Rossum, G., &amp; Drake, F. L. (2009). <i>Python 3 Reference Manual</i>. CreateSpace. &middot;
              <a href="https://www.python.org" target="_blank" rel="noopener">python.org</a></td>
          </tr>
          <tr>
            <td>Nextstrain</td><td class="mono">데이터 출처 (버전 무관, 상시 최신)</td>
            <td>LASV S/L segment population 서열 + lineage 메타데이터 + 계통수 원본</td>
            <td>Hadfield, J., et al. (2018). Nextstrain: real-time tracking of pathogen evolution.
              <i>Bioinformatics</i>, 34(23), 4121&ndash;4123.
              <a href="https://doi.org/10.1093/bioinformatics/bty407" target="_blank" rel="noopener">doi:10.1093/bioinformatics/bty407</a></td>
          </tr>
        </tbody>
      </table>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>분석 파라미터</h3>
      <div class="guide-params">
        <div class="guide-param-row"><span class="guide-param-name">Primer mismatch 허용치 (기본 1)</span><span class="guide-param-desc">F/R 프라이머가 게놈에 결합했다고 인정할 최대 mismatch 개수. mfeprimer는 0/1/2 세 단계를 모두 미리 계산해두고, 이 값으로 판정 기준(reconcile)만 정합니다.</span></div>
        <div class="guide-param-row"><span class="guide-param-name">Probe mismatch 허용치 (기본 0, 완전일치)</span><span class="guide-param-desc">TaqMan probe는 원리상 mismatch에 특히 민감하므로 기본값은 완전일치입니다.</span></div>
        <div class="guide-param-row"><span class="guide-param-name">목표 Coverage (%, 기본 95.0)</span><span class="guide-param-desc">이 값 이상이면 "충분히 커버함"으로 표시하는 참고 기준선입니다.</span></div>
        <div class="guide-param-row"><span class="guide-param-name">3' 말단 연속 매치 요구 염기수 (기본 2)</span><span class="guide-param-desc">프라이머 3' 말단 쪽 mismatch가 폴리머라제 신장 효율에 특히 치명적이라는 점을 반영한 값 — F/R 각각의 3' 말단 N염기 중 하나라도 mismatch가 있으면, 전체 mismatch 개수가 허용치 이내라도 검증 실패로 처리합니다. 값을 바꾸면 mfeprimer 결과를 다시 계산해야 해서 "분석 다시하기"가 필요합니다.</span></div>
      </div>
      <div class="note" style="margin-top:8px;">"⚙ 분석 파라미터 설정" 페이지에서 값을 바꾸고 "분석 다시하기"를 누르면 반영됩니다.</div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>결과 해석 가이드</h3>
      <div class="guide-legend-grid">
        <div class="guide-legend-item"><span class="msa-mm-tag mm-0">0</span> mismatch 없음 (완전일치)</div>
        <div class="guide-legend-item"><span class="msa-mm-tag mm-1">1</span> mismatch 1개</div>
        <div class="guide-legend-item"><span class="msa-mm-tag mm-2plus">2+</span> mismatch 2개 이상</div>
        <div class="guide-legend-item"><span class="pill pill-ok">O 매치</span> / <span class="pill pill-bad">X 불일치</span> 표준물질(KDCA) 매치 여부</div>
        <div class="guide-legend-item">Tm (결합 온도, &deg;C) — 높을수록 결합이 안정적. 일반적으로 primer 55~65&deg;C, probe 63~72&deg;C 권장.</div>
        <div class="guide-legend-item">ΔG (결합 자유에너지, kcal/mol) — 더 음수일수록(예: -20 이하) 결합이 안정적. 0에 가까울수록 결합이 약함.</div>
        <div class="guide-legend-item">amplicon 크기 — TaqMan qPCR은 보통 60~150bp가 이상적이며, 클수록 증폭 효율이 떨어질 수 있음.</div>
        <div class="guide-legend-item">⚠ 교차 결합 — 그룹 내 서로 다른 세트의 F/R/Probe가 의도치 않게 결합하는 경우.</div>
      </div>
      <div class="note" style="margin-top:8px;">각 게놈/표준물질/교차 결합 상세 보기에는 이 기준을 바탕으로 한 "검증 통과/실패" 요약과 "결합 품질 평가(참고용)"가 자동으로 표시됩니다.</div>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>ΔG (결합 자유에너지)란? 어떻게 해석하나요?</h3>
      <p>ΔG(Gibbs free energy of hybridization)는 프라이머/프로브가 상보적인 서열과 이중나선을 이룰 때
      <b>그 결합이 얼마나 자발적이고 안정적으로 일어나는지</b>를 나타내는 열역학적 지표입니다(단위:
      kcal/mol). mfeprimer가 각 F/Probe/R 결합마다 실제 서열 조성(길이, GC 함량, 어떤 염기끼리
      마주보는지)을 바탕으로 계산합니다.</p>
      <ul class="guide-stage-items" style="margin-top:8px;">
        <li><b>음수가 클수록(더 음수, 예: -25 kcal/mol) 결합이 강하고 안정적</b>입니다 — 에너지를
        방출하며 자발적으로 결합이 일어난다는 뜻입니다.</li>
        <li><b>0에 가까울수록(약한 음수, 예: -8 kcal/mol) 결합이 약합니다</b> — annealing 온도에서
        쉽게 떨어지거나(해리), 애초에 잘 결합하지 않을 수 있습니다.</li>
        <li>이 대시보드에서는 대략 <b>-14 kcal/mol보다 약한(0에 더 가까운) 경우</b>를 "결합 품질 평가"에서
        ⚠️로 표시합니다 — 절대적인 기준이라기보다, 같은 세트 안에서 유독 약한 결합을 찾아내기 위한
        참고선입니다.</li>
        <li><b>Tm과의 관계</b>: ΔG와 Tm은 같은 결합 안정성을 다른 방식으로 보여줍니다 — Tm이 "몇 도에서
        풀리는가"라면, ΔG는 "그 온도(보통 mfeprimer는 annealing 온도 근처를 기준으로 계산)에서 결합이
        얼마나 우세한가"입니다. 보통 Tm이 높으면 ΔG도 더 음수인 경향이 있지만, 항상 정비례하지는
        않습니다.</li>
        <li><b>mismatch=0인데 ΔG가 약하다면</b>: 완전일치라도 프라이머가 짧거나 AT(A/T) 비율이 높으면
        ΔG가 약해질 수 있습니다 — "서열은 정확히 일치하지만 결합력 자체가 약해서 실제 PCR에서
        비효율적일 수 있다"는 뜻이므로, mismatch 개수만 보고 안심하기보다 ΔG/Tm도 함께 확인하는 것이
        좋습니다.</li>
      </ul>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>리니지별 유전체 다양성(pixy &#960;)이란? 어떻게 해석하나요?</h3>
      <p>&#960;(파이, nucleotide diversity)는 <b>같은 리니지에 속한 게놈들을 둘씩 짝지어 비교했을 때,
      한 위치(또는 window)에서 평균적으로 염기가 다를 확률</b>입니다. 0이면 그 리니지 안에서는 해당
      위치가 완전히 보존되어 있다는 뜻이고, 값이 클수록(이론상 최대 1에 가까울수록) 그 리니지 내에서도
      변이가 많이 일어난 자리라는 뜻입니다. 이 도구는 population 전체가 아니라 <b>리니지별로 따로</b>
      &#960;를 계산합니다 — "LASV 전체에서 다양하다"와 "내가 검증하려는 특정 리니지 안에서도 다양하다"는
      전혀 다른 의미이기 때문입니다.</p>
      <ul class="guide-stage-items" style="margin-top:8px;">
        <li><b>왜 population 전체가 아니라 mfeprimer 결과처럼 리니지를 나눠서 볼까요?</b> 프라이머가
        "population의 91%를 커버한다"고 해도, 그 91%가 여러 리니지에 고르게 퍼져 있는지 아니면 특정
        리니지에서만 잘 맞고 다른 리니지는 애초에 그 위치가 원래 다양해서(&#960;가 높아서) 놓치기 쉬운
        구조인지는 커버리지 숫자만으로는 알 수 없습니다. 리니지별 &#960;는 "이 결합 위치가 특정 리니지
        안에서도 보존된 자리인가"를 직접 보여줍니다.</li>
        <li><b>전체 LASV(ALL_LASV) 점선 기준선</b>: 리니지별 곡선과 별개로, 필터링된 population 게놈
        전체(리니지 구분 없이 모두 하나의 population으로 취급)를 대상으로 pixy를 한 번 더 돌려 만든
        기준선입니다. 이것은 리니지별 &#960; 값들을 단순 평균한 것과 다릅니다 — 리니지 내부 비교만 반영하는
        개별 곡선과 달리, 이 점선은 <b>리니지 사이의 차이(inter-lineage divergence)</b>까지 포함하기
        때문에 대체로 개별 리니지보다 높거나 중간값 부근에 위치합니다. "이 위치가 LASV 전체를 통틀어
        봤을 때 보존된 자리인가"를 한 눈에 보려면 이 점선을, 특정 리니지 안에서의 보존 여부가 궁금하면
        색깔 있는 실선을 보면 됩니다.</li>
        <li><b>계산 방법 (whole-genome 정렬 → VCF → pixy)</b>: mfeprimer 매치는 프라이머 길이(~20bp)
        단위의 짧은 국소 정렬이라 population 전체의 다양성을 계산하기엔 신호가 너무 적습니다. 그래서
        이 계산만은 별도 경로를 씁니다 — S/L segment의 population 게놈 전체(수백 개)를 <b>mafft
        --adjustdirection</b>로 한 번에 정렬하고(짧은 프라이머 정렬에서는 이 옵션이 신뢰도가 낮아 쓰지
        않지만, 유전체 전체 길이에서는 안정적으로 동작함을 확인했습니다), 참조유전체(S=AY628203,
        L=AY628202) 자신의 정렬된 행을 기준 좌표로 삼아 all-sites VCF(변이 위치뿐 아니라 모든 위치를
        기록)를 직접 구성한 뒤, <a href="https://github.com/ksamuk/pixy" target="_blank" rel="noopener">pixy</a>의
        <code>--stats pi</code>로 리니지별·윈도우별 &#960;를 계산합니다.</li>
        <li><b>왜 굳이 --adjustdirection이 필요한가요?</b> GenBank에 등록된 LASV S/L 게놈 중 일부는
        반대 방향(strand)으로 제출되어 있습니다 — 단순히 참조유전체에 그대로 정렬하면 같은 리니지인데도
        절반 가까이가 정렬 자체에 실패하거나 비정상적으로 높은 mismatch(예: 약 55%)로 나타납니다. 방향을
        맞춰주면 실제 LASV 리니지 간 다양성 수준(약 20% 안팎)과 맞는 정상적인 정렬 결과가 나옵니다.</li>
        <li><b>윈도우 크기 1bp(염기 단위)인 이유</b>: pixy는 보통 몇십&ndash;몇백 bp 단위 window로
        genome-wide &#960;를 요약하지만(Nextstrain의 genome-wide entropy 패널도 비슷한 개념), 이
        도구가 실제로 보고 싶은 것은 <b>하나의 amplicon(보통 100&ndash;170bp) 안에서, 심지어 한 염기
        단위로</b> 결합 위치별 다양성이 어떻게 분포하는지입니다. 50bp 단위로는 amplicon 하나에 막대가
        2~3개만 나와 구조를 알아보기 어려워서 염기 단위(1bp)까지 촘촘하게 잡았습니다. 그 대신 window가
        1bp면 대부분의 위치는 값이 없거나(no_sites=0, 그 리니지 표본들이 전부 결측인 위치) 0(완전
        보존)이고, 실제 변이가 있는 정확한 염기에서만 뾰족하게 튀는 "스파이크" 형태가 됩니다 — 리니지별
        오버레이 차트는 이 스파이크 형태 그대로 보여주고, MSA 위 히스토그램은 여러 리니지를 평균 내기
        때문에 상대적으로 덜 들쭉날쭉합니다. 즉 "이 구간이 대체로 얼마나 다양한가"보다 "정확히 이 염기가
        어떤 리니지에서 변이하는가"를 짚어내는 데 더 적합한 해상도입니다.</li>
        <li><b>어디서 볼 수 있나요?</b> "🧬 리니지별 유전체 다양성" 페이지에서는 segment 전체 길이에 대해
        리니지별 &#960; 곡선을 유전자 경계와 함께 볼 수 있고, 각 프라이머 세트 상세 페이지에서는 그 세트의
        F/Probe/R 결합 위치 주변으로 확대한 오버레이 차트, 그리고 전체 게놈 정렬(MSA) 화면의 position
        룰러 바로 아래에 같은 데이터를 막대그래프(히스토그램)로 표시합니다.</li>
        <li><b>MSA 위 히스토그램의 근사치 주의사항</b>: 전체 게놈 정렬(MSA)은 이 세트만의 로컬 정렬이라
        pixy의 참조유전체 좌표계와 완전히 1:1로 맞지는 않습니다 — 세트의 참조유전체 결합위치(예:
        AY628203:1222)를 기준점 삼아 컬럼 번호를 근사적으로 매핑한 것이라, 실제 위치와 몇 bp 정도
        어긋날 수 있습니다(게놈 대부분이 짧은 indel 없이 정렬되는 amplicon 구간에서는 오차가 거의
        없습니다). 정확한 좌표가 필요하면 세트 상세 페이지의 오버레이 차트나 "🧬 리니지별 유전체 다양성"
        페이지를 참고하세요.</li>
        <li><b>이 값이 낮다고 프라이머가 무조건 좋은 것은 아닙니다</b>: &#960;는 population에 이미
        존재하는 다양성만 반영합니다 — 표본으로 확보된 게놈 수가 적은 리니지(예: LI, LVII)는 실제보다
        &#960;가 낮게 나올 수 있으므로, 표준물질 매치·게놈별 mismatch 표 등 다른 지표와 함께 보는 것이
        좋습니다.</li>
      </ul>
    </div>

    <div class="card" style="margin-bottom:20px;">
      <h3>용어 설명</h3>
      ${renderGuideFootnoteList()}
    </div>
  `;
}

function showGuidePage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("guideView");
  $("#guideView").innerHTML = renderGuidePage();
}

$("#guideBtn").addEventListener("click", showGuidePage);

// ---------------------------------------------------------------------------------------
// 리니지별 유전체 다양성 (pixy pi) -- Nextstrain의 genome-wide entropy 패널과 유사하게,
// S/L segment 각각에 대해 리니지별 pi를 genome position에 따른 line chart로 표시한다.
const DIVERSITY_PALETTE = ["#4C78A8", "#E45756", "#54A24B", "#F2A65A", "#B279A2", "#72B7B2", "#FF9DA6", "#9D755D", "#BAB0AC", "#EECA3B"];
function colorForIndex(i) { return DIVERSITY_PALETTE[i % DIVERSITY_PALETTE.length]; }

// run_pixy.py also computes one extra "population" that pools every sample regardless of
// lineage -- not the same as averaging the per-lineage curves (which only counts
// within-lineage pairs), since pooling also captures inter-lineage divergence. Rendered as a
// distinct dashed baseline rather than just another colored line in the lineage palette.
// The virus profile (profiles/<id>.json via /api/profile) names the pooled all-genomes population and the
// short virus label; these defaults are replaced at start-up by loadProfile().
let PROFILE = { id: "lassa", name: "Lassa mammarenavirus", short_name: "LASV", population: { pooled_label: "ALL_LASV" } };
let ALL_LASV_KEY = PROFILE.population.pooled_label;
async function loadProfile() {
  try {
    PROFILE = await fetchJSON("/api/profile");
    ALL_LASV_KEY = PROFILE.population.pooled_label;
    document.title = `${PROFILE.short_name} Primer/Probe Validation`;
    const fv = document.getElementById("footerVersion");
    if (fv && PROFILE.version) fv.textContent = PROFILE.version;
    const h1 = document.querySelector(".sidebar-header h1");
    if (h1) h1.textContent = `${PROFILE.short_name} Primer/Probe`;
  } catch { /* keep the defaults */ }
  renderTaxonLine();
  initPipelineSwitch();
  try {                                                   // the figures built in the browser must follow the same pass rule as the server's
    const prm = await fetchJSON("/api/params");
    if (Number.isFinite(prm.primer_mismatch_pass)) PASS.mm = prm.primer_mismatch_pass;
    if (Number.isFinite(prm.probe_mismatch)) PASS.pm = prm.probe_mismatch;
  } catch { /* keep the defaults */ }
}

// ---- ICTV taxonomy (profiles/<id>.json "taxonomy"): species under the title, the full rank list in a pop-up ----
function ictvCrumbHtml() {
  const t = PROFILE.taxonomy;
  if (!t || !t.ranks) return "";
  return `<div class="ictv-crumb" title="${escHtml(t.source || "")}">ICTV 분류: ${t.ranks.map(([rk, nm], i) => (i === t.ranks.length - 1 ? `<i><b>${escHtml(nm)}</b></i>` : `<a class="ictv-open" title="${escHtml(rk)}">${escHtml(nm)}</a>`) + (i === t.ranks.length - 1 ? ` <span>(${escHtml(rk)})</span>` : "")).join(" &rsaquo; ")}</div>`;
}
function renderTaxonLine() {
  const b = document.getElementById("taxonLine"), t = PROFILE.taxonomy;
  if (!b || !t || !t.ranks) return;
  const species = t.ranks[t.ranks.length - 1][1], genus = (t.ranks.find((r) => r[0] === "Genus") || [])[1];
  b.innerHTML = `ICTV · <i>${escHtml(species)}</i> <span style="opacity:.7;">(${escHtml(genus || "")})</span> &#9656;`;
  b.hidden = false;
  b.onclick = openTaxonomyModal;
}
function openTaxonomyModal() {
  const t = PROFILE.taxonomy;
  if (!t) return;
  openDetailModal({
    title: `ICTV 분류 &mdash; ${escHtml(t.virus || PROFILE.short_name)}`, subtitle: "Realm &rsaquo; &hellip; &rsaquo; Family &rsaquo; Genus &rsaquo; Species", compact: true,
    html: `<div class="ictv-ranks">${t.ranks.map(([rk, nm], i) => `<div class="rk">${escHtml(rk)}</div><div class="nm" style="padding-left:${Math.min(i, 9) * 6}px;">${i === t.ranks.length - 1 ? `<i>${escHtml(nm)}</i>` : escHtml(nm)}</div>`).join("")}</div>
      ${t.former_species ? `<p class="note">이전 종 이름: <i>${escHtml(t.former_species)}</i></p>` : ""}
      <p class="note">${escHtml(t.source || "")}</p>`,
  });
}
// the LASV / NiV switch in the sidebar: the dashboards listed in entry/pipelines.json, opened on the same host
async function initPipelineSwitch() {
  const box = document.getElementById("pipelineSwitch");
  const link = document.getElementById("entryLink");
  if (link) link.href = `${location.protocol}//${location.hostname}:8760/`;
  if (!box) return;
  let list = [];
  try { list = await fetchJSON("/api/pipelines"); } catch (e) { /* stand-alone project */ }
  if (list.length < 2) return;
  box.innerHTML = list.map((p) => (p.current
    ? `<span class="pipe-pill current" title="지금 보고 있는 대시보드">${escHtml(p.name)}</span>`
    : `<a class="pipe-pill" href="${location.protocol}//${location.hostname}:${p.port}/" title="${escHtml(p.name)} 대시보드로 이동">${escHtml(p.name)}</a>`)).join("");
  box.hidden = false;
}
document.addEventListener("click", (e) => { if (e.target.closest(".ictv-open")) openTaxonomyModal(); });
function splitLineageKeys(segData) {
  const keys = Object.keys(segData.windows || {}).sort();
  return { real: keys.filter((k) => k !== ALL_LASV_KEY), hasAll: keys.includes(ALL_LASV_KEY) };
}

const diversityState = { seg: null, data: null, hidden: new Set(), xDomain: null, smooth: 25 };

// pixy at window=1bp gives one value per site, which is far too spiky to read with 5-8 lineages
// overlaid -- so the main diversity page smooths each curve with a centered moving average over
// `smooth` bp. It is weighted by pixy's own per-window no_sites (sum(pi*no_sites)/sum(no_sites)),
// i.e. the same weighting pixy uses for a genome-wide figure, so a stretch with few comparable
// sites doesn't count as much as a well-covered one. smooth=1 shows the raw per-site values.
const _smoothCache = new WeakMap();
function smoothDiversitySeries(windows, smoothBp, windowSize) {
  if (!windows || smoothBp <= 1) return windows || [];
  const k = Math.max(1, Math.round(smoothBp / Math.max(1, windowSize || 1)));
  if (k <= 1) return windows;
  let byK = _smoothCache.get(windows);
  if (!byK) { byK = new Map(); _smoothCache.set(windows, byK); }
  if (byK.has(k)) return byK.get(k);
  const n = windows.length;
  const num = new Float64Array(n + 1), den = new Float64Array(n + 1);
  for (let i = 0; i < n; i++) {
    const w = windows[i];
    const ok = w.pi != null && w.no_sites > 0;
    num[i + 1] = num[i] + (ok ? w.pi * w.no_sites : 0);
    den[i + 1] = den[i] + (ok ? w.no_sites : 0);
  }
  const h = Math.floor(k / 2);
  const out = new Array(n);
  for (let i = 0; i < n; i++) {
    const lo = Math.max(0, i - h), hi = Math.min(n, i + h + 1);
    const d = den[hi] - den[lo];
    out[i] = { start: windows[i].start, end: windows[i].end, no_sites: windows[i].no_sites,
               pi: d > 0 ? (num[hi] - num[lo]) / d : null };
  }
  byK.set(k, out);
  return out;
}

function renderDiversityGeneTrack(genes, xScale, trackY, trackH, plotLeft, plotRight) {
  if (!genes || !genes.length) return "";
  return genes.map((g) => {
    const x0 = Math.max(plotLeft, xScale(g.start));
    const x1 = Math.min(plotRight, xScale(g.end));
    if (x1 <= x0) return "";
    const w = x1 - x0;
    const fill = g.strand === "-" ? "#9aa3af" : "#6b7280";
    const label = w > 22 ? `<text x="${(x0 + x1) / 2}" y="${trackY + trackH / 2 + 4}" text-anchor="middle" font-size="10" fill="#fff" font-weight="600">${escHtml(g.gene)}</text>` : "";
    return `<g>
      <rect x="${x0}" y="${trackY}" width="${w}" height="${trackH}" rx="3" fill="${fill}" opacity="0.85"></rect>
      ${label}
    </g>`;
  }).join("");
}

function renderDiversityOverlays(overlays, xScale, plotTop, plotH, plotLeft, plotRight) {
  if (!overlays || !overlays.length) return "";
  // Every label used to sit on the exact same row (plotTop+12) -- fine for one set's 3
  // widely-spaced F/P/R boxes, but a GROUP's combined overlay can have a dozen boxes from
  // different member sets clustered in the same amplicon region, so their labels piled up
  // unreadably on top of each other. Sort by x position and greedily assign each label to
  // the first row (of a small fixed set) whose last-placed label doesn't horizontally
  // overlap it -- a simple, index-free way to de-collide text without real DOM text-width
  // measurement (this is a hand-rolled SVG string template, not a live-measurable canvas).
  const ROW_HEIGHT = 10;
  const MAX_ROWS = 8;
  const items = overlays.map((o) => {
    const x0 = Math.max(plotLeft, xScale(Math.min(o.start, o.end)));
    const x1 = Math.min(plotRight, xScale(Math.max(o.start, o.end)));
    const labelX = Math.min(Math.max((x0 + x1) / 2, plotLeft + 20), plotRight - 20);
    // Rough label half-width estimate (px) from character count, since we can't measure
    // rendered text width here -- generous enough to avoid most false "no overlap" calls.
    const halfWidth = Math.max(10, o.label.length * 3.4);
    return { o, x0, x1, labelX, labelLeft: labelX - halfWidth, labelRight: labelX + halfWidth };
  }).filter((it) => it.x1 > it.x0);
  items.sort((a, b) => a.labelLeft - b.labelLeft);
  const rowLastRight = new Array(MAX_ROWS).fill(-Infinity);
  items.forEach((it) => {
    let row = rowLastRight.findIndex((right) => it.labelLeft > right + 4);
    if (row === -1) row = items.indexOf(it) % MAX_ROWS; // every row already collides -- just cycle
    rowLastRight[row] = it.labelRight;
    it.row = row;
  });
  return items.map(({ o, x0, x1, labelX, row }) => {
    const w = x1 - x0;
    return `<g>
      <rect x="${x0.toFixed(1)}" y="${plotTop}" width="${w.toFixed(1)}" height="${plotH}" fill="${o.color}" opacity="0.15"></rect>
      <line x1="${x0.toFixed(1)}" y1="${plotTop}" x2="${x0.toFixed(1)}" y2="${plotTop + plotH}" stroke="${o.color}" stroke-width="1" stroke-dasharray="3,2" opacity="0.7"></line>
      <line x1="${x1.toFixed(1)}" y1="${plotTop}" x2="${x1.toFixed(1)}" y2="${plotTop + plotH}" stroke="${o.color}" stroke-width="1" stroke-dasharray="3,2" opacity="0.7"></line>
      <text x="${labelX.toFixed(1)}" y="${plotTop + 12 + row * ROW_HEIGHT}" text-anchor="middle" font-size="10" font-weight="700" fill="${o.color}">${escHtml(o.label)}</text>
    </g>`;
  }).join("");
}

function renderDiversityChart(segKey, segData, genes, genomeLen, opts) {
  opts = opts || {};
  const xMin = opts.xMin != null ? opts.xMin : 1;
  const xMax = opts.xMax != null ? opts.xMax : genomeLen;
  const overlays = opts.overlays || [];
  const width = 920, height = 300;
  const marginLeft = 54, marginRight = 16, marginTop = 16, marginBottom = 56;
  const geneTrackH = 22;
  const plotW = width - marginLeft - marginRight;
  const plotTop = marginTop, plotBottom = height - marginBottom - geneTrackH - 8;
  const plotH = plotBottom - plotTop;
  const plotLeft = marginLeft, plotRight = width - marginRight;

  const { real: lineages, hasAll } = splitLineageKeys(segData);
  const allKeys = hasAll ? lineages.concat([ALL_LASV_KEY]) : lineages;
  let maxPi = 0.01;
  for (const lin of allKeys) {
    for (const w of segData.windows[lin]) {
      if (w.pi != null && w.start <= xMax && w.end >= xMin) maxPi = Math.max(maxPi, w.pi);
    }
  }
  maxPi *= 1.1;

  const xScale = (pos) => marginLeft + ((pos - xMin) / (xMax - xMin)) * plotW;
  const yScale = (pi) => plotTop + plotH - (pi / maxPi) * plotH;

  const yTicks = 5;
  let gridSvg = "";
  for (let t = 0; t <= yTicks; t++) {
    const pi = (maxPi / yTicks) * t;
    const y = yScale(pi);
    gridSvg += `<line x1="${marginLeft}" y1="${y.toFixed(1)}" x2="${width - marginRight}" y2="${y.toFixed(1)}" stroke="var(--border)" stroke-width="1"></line>`;
    gridSvg += `<text x="${marginLeft - 8}" y="${(y + 3).toFixed(1)}" text-anchor="end" font-size="10" fill="var(--muted)">${pi.toFixed(3)}</text>`;
  }
  const xTicks = 6;
  let xAxisSvg = "";
  for (let t = 0; t <= xTicks; t++) {
    const pos = Math.round(xMin + ((xMax - xMin) / xTicks) * t);
    const x = xScale(pos);
    xAxisSvg += `<line x1="${x.toFixed(1)}" y1="${plotTop}" x2="${x.toFixed(1)}" y2="${plotBottom}" stroke="var(--border)" stroke-width="1" opacity="0.4"></line>`;
    xAxisSvg += `<text x="${x.toFixed(1)}" y="${plotBottom + 16}" text-anchor="middle" font-size="10" fill="var(--muted)">${pos}</text>`;
  }

  let linesSvg = "";
  lineages.forEach((lin, i) => {
    const color = colorForIndex(i);
    const windows = segData.windows[lin];
    let d = "";
    let started = false;
    for (const w of windows) {
      if (w.pi == null || w.end < xMin || w.start > xMax) { started = false; continue; }
      const x = xScale((w.start + w.end) / 2);
      const y = yScale(w.pi);
      d += (started ? " L " : " M ") + x.toFixed(1) + "," + y.toFixed(1);
      started = true;
    }
    linesSvg += `<path d="${d}" fill="none" stroke="${color}" stroke-width="1.6" opacity="0.9"></path>`;
  });
  if (hasAll) {
    const windows = segData.windows[ALL_LASV_KEY];
    let d = "";
    let started = false;
    for (const w of windows) {
      if (w.pi == null || w.end < xMin || w.start > xMax) { started = false; continue; }
      const x = xScale((w.start + w.end) / 2);
      const y = yScale(w.pi);
      d += (started ? " L " : " M ") + x.toFixed(1) + "," + y.toFixed(1);
      started = true;
    }
    linesSvg += `<path d="${d}" fill="none" stroke="var(--fg)" stroke-width="2.2" stroke-dasharray="7,4" opacity="0.75"></path>`;
  }

  const overlaySvg = renderDiversityOverlays(overlays, xScale, plotTop, plotH, plotLeft, plotRight);
  const geneTrackY = plotBottom + 26;
  const geneSvg = renderDiversityGeneTrack(genes, xScale, geneTrackY, geneTrackH, plotLeft, plotRight);

  return `<svg viewBox="0 0 ${width} ${height}" class="diversity-svg" preserveAspectRatio="xMidYMid meet">
    ${gridSvg}
    ${xAxisSvg}
    ${geneSvg}
    ${overlaySvg}
    ${linesSvg}
    <text x="${marginLeft}" y="12" font-size="10" fill="var(--muted)">&#960; (nucleotide diversity, window=${segData.window_size}bp)</text>
    <text x="${width - marginRight}" y="${geneTrackY + geneTrackH + 14}" text-anchor="end" font-size="10" fill="var(--muted)">genome position (nt, ${escHtml(segKey.toUpperCase())} segment, ref length ${genomeLen}bp)</text>
  </svg>`;
}

function renderDiversityLegend(lineages, genomeWidePi, lineageTotals) {
  lineageTotals = lineageTotals || {};
  const hasAll = lineages.includes(ALL_LASV_KEY);
  const realLineages = lineages.filter((l) => l !== ALL_LASV_KEY);
  const allItem = hasAll ? (() => {
    const avg = genomeWidePi[ALL_LASV_KEY];
    return `<div class="diversity-legend-item diversity-legend-all"><span class="swatch swatch-dashed"></span><b>전체 ${PROFILE.short_name} (population 전체, 리니지 혼합)</b>
      <span class="muted">(genome-wide &#960; ${avg != null ? avg.toFixed(4) : "N/A"})</span></div>`;
  })() : "";
  const lineageItems = realLineages.map((lin, i) => {
    const color = colorForIndex(i);
    const avg = genomeWidePi[lin];
    const n = lineageTotals[lin];
    return `<div class="diversity-legend-item"><span class="swatch" style="background:${color}"></span>${escHtml(lin)}
      <span class="muted">(genome-wide &#960; ${avg != null ? avg.toFixed(4) : "N/A"})</span> ${lowNBadge(n)}</div>`;
  }).join("");
  return `<div class="diversity-legend">${allItem}${lineageItems}</div>`;
}

// Interactive (zoomable + toggleable) version of the main "리니지별 유전체 다양성" page's chart,
// built with d3 (brush-to-zoom + a legend viewers can click to hide/show a lineage's line) --
// kept separate from renderDiversityChart() above, which stays a plain string-rendered SVG
// because it's also reused by the per-set binding-site overlay (loadSetDiversityOverlay), a
// fixed, non-interactive zoomed-in window that has no need for drag-to-zoom or a toggle legend.
// diversityState.hidden (a Set of lineage keys) and diversityState.xDomain (the current
// zoomed-in {seg, range}) persist across re-renders so a segment-tab switch or legend click
// doesn't reset what the viewer was just looking at (xDomain resets on its own when the segment
// doesn't match, since S/L have different genome lengths).
function renderInteractiveDiversityChart(container, segKey, segData, genes, genomeLen) {
  container.innerHTML = "";
  const { real: lineages, hasAll } = splitLineageKeys(segData);
  const fullDomain = [1, genomeLen];
  let xDomain = (diversityState.xDomain && diversityState.xDomain.seg === segKey)
    ? diversityState.xDomain.range : fullDomain;

  const toolbar = document.createElement("div");
  toolbar.className = "diversity-toolbar";
  toolbar.innerHTML = `<button id="diversityResetZoomBtn" class="btn btn-secondary">&#8634; 전체 보기</button>
    <label class="diversity-smooth-label">평활화
      <select id="diversitySmoothSel">${[1, 10, 25, 50, 100].map((v) => `<option value="${v}" ${v === diversityState.smooth ? "selected" : ""}>${v === 1 ? "없음 (1bp 원본)" : `${v}bp 이동평균`}</option>`).join("")}</select>
    </label>
    <span class="note">차트를 드래그해 영역을 선택하면 확대됩니다 &middot; 더블클릭하면 초기화됩니다 &middot; 마우스를 올리면 위치별 π를 볼 수 있습니다.</span>
    <div class="diversity-download-row">
      <label>시작 <input id="diversityDlStart" type="number" min="1" max="${genomeLen}"></label>
      <label>끝 <input id="diversityDlEnd" type="number" min="1" max="${genomeLen}"></label>
      <select id="diversityDlLineage"><option value="${ALL_LASV_KEY}">전체 ${PROFILE.short_name}</option>${lineages.map((lin) => `<option value="${escHtml(lin)}">${escHtml(lin)}</option>`).join("")}</select>
      <button id="diversityDlBtn" class="btn btn-secondary" title="지정한 참조유전체 위치 구간에 대해, 필터링된 각 게놈이 실제로 갖는 서열을 FASTA로 받습니다">&#8681; 이 구간 서열 다운로드</button>
      <button id="diversityAlignBtn" class="btn btn-secondary" title="지정한 참조유전체 위치 구간에 대해, 필터링된 각 게놈의 실제 서열을 아래에 정렬(alignment)로 보여줍니다">&#129516; 이 구간 정렬 보기</button>
    </div>`;
  container.appendChild(toolbar);

  const dlStart = toolbar.querySelector("#diversityDlStart");
  const dlEnd = toolbar.querySelector("#diversityDlEnd");
  const syncDlInputs = () => { dlStart.value = xDomain[0]; dlEnd.value = xDomain[1]; };
  syncDlInputs();
  const readRange = () => {
    const lo = parseInt(dlStart.value, 10), hi = parseInt(dlEnd.value, 10);
    if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo < 1 || hi < lo) {
      uiAlert("시작/끝 위치를 올바르게 입력해주세요 (시작 ≤ 끝, 1 이상).");
      return null;
    }
    return { lo, hi, lineage: toolbar.querySelector("#diversityDlLineage").value };
  };
  toolbar.querySelector("#diversityDlBtn").addEventListener("click", () => {
    const r = readRange();
    if (!r) return;
    window.open(`/api/diversity/download?seg=${segKey}&start=${r.lo}&end=${r.hi}&lineage=${encodeURIComponent(r.lineage)}`, "_blank");
  });
  toolbar.querySelector("#diversityAlignBtn").addEventListener("click", () => {
    const r = readRange();
    if (!r) return;
    loadDiversityRegionAlignment(segKey, r.lo, r.hi, r.lineage);
  });

  const width = 920, height = 300;
  const marginLeft = 54, marginRight = 16, marginTop = 16, marginBottom = 56;
  const geneTrackH = 22;
  const plotLeft = marginLeft, plotRight = width - marginRight;
  const plotTop = marginTop, plotBottom = height - marginBottom - geneTrackH - 8;
  const geneTrackY = plotBottom + 26;

  const svg = d3.select(container).append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`).attr("class", "diversity-svg");
  const seriesFor = (lin) => smoothDiversitySeries(segData.windows[lin] || [], diversityState.smooth, segData.window_size);
  toolbar.querySelector("#diversitySmoothSel").addEventListener("change", (ev) => {
    diversityState.smooth = parseInt(ev.target.value, 10) || 1;
    render();
  });
  const tip = document.createElement("div");
  tip.className = "diversity-tip";
  tip.hidden = true;
  container.style.position = "relative";

  function render() {
    svg.selectAll("*").remove();
    const [xMin, xMax] = xDomain;
    const x = d3.scaleLinear().domain([xMin, xMax]).range([plotLeft, plotRight]);

    const visibleKeys = lineages.concat(hasAll ? [ALL_LASV_KEY] : []).filter((k) => !diversityState.hidden.has(k));
    let maxPi = 0.01;
    for (const lin of visibleKeys) {
      for (const w of seriesFor(lin)) {
        if (w.pi != null && w.start <= xMax && w.end >= xMin) maxPi = Math.max(maxPi, w.pi);
      }
    }
    maxPi *= 1.1;
    const y = d3.scaleLinear().domain([0, maxPi]).range([plotBottom, plotTop]);

    const yTicks = 5;
    for (let t = 0; t <= yTicks; t++) {
      const pi = (maxPi / yTicks) * t;
      const yy = y(pi);
      svg.append("line").attr("x1", plotLeft).attr("y1", yy).attr("x2", plotRight).attr("y2", yy).attr("stroke", "var(--border)");
      svg.append("text").attr("x", plotLeft - 8).attr("y", yy + 3).attr("text-anchor", "end").attr("font-size", 10).attr("fill", "var(--muted)").text(pi.toFixed(3));
    }
    const xTicks = 6;
    for (let t = 0; t <= xTicks; t++) {
      const pos = Math.round(xMin + ((xMax - xMin) / xTicks) * t);
      const xx = x(pos);
      svg.append("line").attr("x1", xx).attr("y1", plotTop).attr("x2", xx).attr("y2", plotBottom).attr("stroke", "var(--border)").attr("opacity", 0.4);
      svg.append("text").attr("x", xx).attr("y", plotBottom + 16).attr("text-anchor", "middle").attr("font-size", 10).attr("fill", "var(--muted)").text(pos);
    }

    (genes || []).forEach((gn) => {
      const x0 = Math.max(plotLeft, x(gn.start)), x1 = Math.min(plotRight, x(gn.end));
      if (x1 <= x0) return;
      svg.append("rect").attr("x", x0).attr("y", geneTrackY).attr("width", x1 - x0).attr("height", geneTrackH).attr("rx", 3)
        .attr("fill", gn.strand === "-" ? "#9aa3af" : "#6b7280").attr("opacity", 0.85);
      if (x1 - x0 > 22) {
        svg.append("text").attr("x", (x0 + x1) / 2).attr("y", geneTrackY + geneTrackH / 2 + 4).attr("text-anchor", "middle")
          .attr("font-size", 10).attr("fill", "#fff").attr("font-weight", 600).text(gn.gene);
      }
    });

    const lineGen = d3.line()
      .defined((d) => d.pi != null && d.end >= xMin && d.start <= xMax)
      .x((d) => x((d.start + d.end) / 2)).y((d) => y(d.pi));
    lineages.forEach((lin, i) => {
      if (diversityState.hidden.has(lin)) return;
      svg.append("path").datum(seriesFor(lin)).attr("fill", "none").attr("stroke", colorForIndex(i))
        .attr("stroke-width", 1.6).attr("opacity", 0.9).attr("d", lineGen);
    });
    if (hasAll && !diversityState.hidden.has(ALL_LASV_KEY)) {
      svg.append("path").datum(seriesFor(ALL_LASV_KEY)).attr("fill", "none").attr("stroke", "var(--fg)")
        .attr("stroke-width", 2.2).attr("stroke-dasharray", "7,4").attr("opacity", 0.75).attr("d", lineGen);
    }

    svg.append("text").attr("x", plotLeft).attr("y", 12).attr("font-size", 10).attr("fill", "var(--muted)")
      .text(diversityState.smooth > 1
        ? `π (nucleotide diversity, ${diversityState.smooth}bp 이동평균)`
        : `π (nucleotide diversity, ${segData.window_size}bp 원본)`);
    svg.append("text").attr("x", plotRight).attr("y", geneTrackY + geneTrackH + 14).attr("text-anchor", "end")
      .attr("font-size", 10).attr("fill", "var(--muted)")
      .text(`genome position (nt, ${segKey.toUpperCase()} segment, ref length ${genomeLen}bp)`);

    // Hover tooltip: nearest position -> every visible curve's smoothed pi there, highest first.
    const guide = svg.append("line").attr("y1", plotTop).attr("y2", plotBottom).attr("stroke", "var(--muted)")
      .attr("stroke-dasharray", "3,3").attr("opacity", 0).attr("pointer-events", "none");
    const valueAt = (series, pos) => {
      if (!series.length) return null;
      let lo = 0, hi = series.length - 1;
      while (lo < hi) { const mid = (lo + hi) >> 1; if (series[mid].end < pos) lo = mid + 1; else hi = mid; }
      const w = series[lo];
      return w && w.start <= pos && w.end >= pos ? w.pi : null;
    };
    const hideTip = () => { tip.hidden = true; guide.attr("opacity", 0); };
    svg.on("mousemove", (event) => {
      const [mx] = d3.pointer(event, svg.node());
      if (mx < plotLeft || mx > plotRight) return hideTip();
      const pos = Math.round(x.invert(mx));
      const rows = [];
      lineages.forEach((lin, i) => {
        if (diversityState.hidden.has(lin)) return;
        const v = valueAt(seriesFor(lin), pos);
        if (v != null) rows.push({ label: lin, color: colorForIndex(i), v, n: (segData.lineage_n || {})[lin] });
      });
      if (hasAll && !diversityState.hidden.has(ALL_LASV_KEY)) {
        const v = valueAt(seriesFor(ALL_LASV_KEY), pos);
        if (v != null) rows.push({ label: `전체 ${PROFILE.short_name}`, color: "var(--fg)", v, n: (segData.lineage_n || {})[ALL_LASV_KEY] });
      }
      rows.sort((a, b) => b.v - a.v);
      const gene = (genes || []).find((gn) => pos >= gn.start && pos <= gn.end);
      tip.innerHTML = `<div class="tip-head">${segKey.toUpperCase()}:${pos}${gene ? ` &middot; ${escHtml(gene.gene)}` : ""}</div>` +
        rows.map((r) => `<div class="tip-row"><span class="dot" style="background:${r.color}"></span>${escHtml(r.label)}${r.n != null ? ` <span class="muted">(n=${r.n})</span>` : ""}<b>${r.v.toFixed(4)}</b></div>`).join("");
      tip.hidden = false;
      const box = container.getBoundingClientRect();
      const px = event.clientX - box.left, py = event.clientY - box.top;
      tip.style.left = `${px + 14 + 200 > box.width ? px - 214 : px + 14}px`;
      tip.style.top = `${Math.max(0, py - 10)}px`;
      guide.attr("x1", mx).attr("x2", mx).attr("opacity", 0.6);
    }).on("mouseleave", hideTip);

    const brush = d3.brushX().extent([[plotLeft, plotTop], [plotRight, plotBottom]])
      .on("end", (event) => {
        if (!event.selection) return;
        const [sx0, sx1] = event.selection.map(x.invert);
        if (sx1 - sx0 < 5) return; // ignore accidental clicks/tiny drags
        xDomain = [Math.max(1, Math.round(sx0)), Math.min(genomeLen, Math.round(sx1))];
        diversityState.xDomain = { seg: segKey, range: xDomain };
        syncDlInputs();
        render();
      });
    svg.append("g").attr("class", "diversity-brush").call(brush)
      .on("dblclick.zoom", () => {
        xDomain = fullDomain;
        diversityState.xDomain = { seg: segKey, range: xDomain };
        syncDlInputs();
        render();
      });
  }

  container.appendChild(tip);
  render();
  container.querySelector("#diversityResetZoomBtn").addEventListener("click", () => {
    xDomain = fullDomain;
    diversityState.xDomain = { seg: segKey, range: xDomain };
    syncDlInputs();
    render();
  });
}

function renderDiversityLegendToggle(lineages, hasAll, genomeWidePi, lineageTotals, lineageN) {
  lineageTotals = lineageTotals || {};
  lineageN = lineageN || {};
  const items = [];
  if (hasAll) {
    const avg = genomeWidePi[ALL_LASV_KEY];
    items.push({ key: ALL_LASV_KEY, color: "var(--fg)",
      label: `전체 ${PROFILE.short_name} (genome-wide &#960; ${avg != null ? avg.toFixed(4) : "N/A"})${lineageN[ALL_LASV_KEY] != null ? ` &middot; n=${lineageN[ALL_LASV_KEY]}` : ""}` });
  }
  lineages.forEach((lin, i) => {
    const avg = genomeWidePi[lin];
    const n = lineageN[lin] != null ? lineageN[lin] : lineageTotals[lin];
    items.push({ key: lin, color: colorForIndex(i),
      label: `${escHtml(lin)} (genome-wide &#960; ${avg != null ? avg.toFixed(4) : "N/A"}) &middot; n=${n != null ? n : "?"}` });
  });
  return `<div class="meta-year-legend">${items.map((it) => {
    const off = diversityState.hidden.has(it.key);
    return `<span class="meta-year-legend-chip${off ? " off" : ""}" data-key="${escHtml(it.key)}" style="--chip-color:${it.color}"><span class="dot"></span>${it.label}</span>`;
  }).join("")}</div>`;
}

// Region-based alignment viewer for the 리니지별 유전체 다양성 page -- reuses the exact same
// .msa-outer/.msa-row/.msa-seq/.msa-label markup as a primer set's own detail-page MSA
// (renderMSA() above) so a manually-picked reference-position window "looks like a detail
// page" even though it isn't tied to any primer set: the data comes from
// /api/diversity/region (server.py's build_region_alignment()), which reconstructs each
// filtered population genome's actual called base at every position straight out of the
// pixy all-sites VCF, the same source the "이 구간 서열 다운로드" FASTA button already uses.
function renderDiversityRegionAlignment(regionData) {
  if (regionData.error) {
    return `<div class="card" style="margin-top:20px;"><p class="note">정렬을 불러올 수 없습니다: ${escHtml(regionData.error)}</p></div>`;
  }
  const { seg, start, end, columns, consensus, rows, total_count, truncated } = regionData;
  const ruler = Array.from({ length: columns }, (_, i) => ((start + i) % 10 === 0 ? "|" : " ")).join("");
  const rowsHtml = rows.map((r) => `<div class="msa-row msa-genome-row">
      <div class="msa-label"><span class="mono">${escHtml(r.genome_id)} | ${escHtml(r.lineage)} | ${escHtml(r.country || "?")}</span></div>
      <div class="msa-seq">${escHtml(r.dots)}</div>
    </div>`).join("");
  return `
    <div class="card" style="margin-top:20px;">
      <h3>${escHtml(seg.toUpperCase())} segment ${start}&ndash;${end} 구간 정렬 (필터링된 게놈 ${total_count}개${truncated ? `, 상위 ${rows.length}개만 표시` : ""})</h3>
      <div class="msa-outer">
        <div class="msa-row msa-ruler"><div class="msa-label">position</div><div class="msa-seq">${ruler}</div></div>
        <div class="msa-row msa-consensus"><div class="msa-label">CONSENSUS</div><div class="msa-seq">${escHtml(consensus)}</div></div>
        ${rowsHtml}
      </div>
      <div class="legend" style="margin-top:8px;">
        <span>"." = population 다수(consensus)와 동일, 글자 = 다수 염기와 다름, "N" = 결측/모호(gap 포함)</span>
        <span>게놈 라벨: ACCESSION | LINEAGE | COUNTRY</span>
      </div>
    </div>`;
}

async function loadDiversityRegionAlignment(seg, lo, hi, lineage) {
  const wrap = document.querySelector("#diversityRegionAlignWrap");
  if (!wrap) return;
  wrap.innerHTML = "<div class='card' style='margin-top:20px;'><p class='note'>정렬 불러오는 중...</p></div>";
  try {
    const res = await fetch(`/api/diversity/region?seg=${seg}&start=${lo}&end=${hi}&lineage=${encodeURIComponent(lineage)}`);
    const data = await res.json();
    wrap.innerHTML = renderDiversityRegionAlignment(data);
  } catch (e) {
    wrap.innerHTML = `<div class="card" style="margin-top:20px;"><p class="note">불러오기 실패: ${escHtml(e.message)}</p></div>`;
  }
}

function mountDiversityChart(view) {
  const chartWrap = view.querySelector("#diversityChartWrap");
  if (!chartWrap) return;
  const data = diversityState.data;
  const seg = diversityState.seg;
  const segData = data.segments[seg];
  const genes = (data.genes || {})[seg] || [];
  const genomeLen = (data.genome_len || {})[seg];
  const { real: lineages, hasAll } = splitLineageKeys(segData);

  renderInteractiveDiversityChart(chartWrap, seg, segData, genes, genomeLen);

  const legendWrap = view.querySelector("#diversityLegendWrap");
  const refreshLegend = () => {
    legendWrap.innerHTML = renderDiversityLegendToggle(lineages, hasAll, segData.genome_wide_pi || {}, (data.lineage_totals || {})[seg] || {}, segData.lineage_n || {});
    legendWrap.querySelectorAll(".meta-year-legend-chip").forEach((chip) => {
      chip.addEventListener("click", () => {
        const key = chip.dataset.key;
        if (diversityState.hidden.has(key)) diversityState.hidden.delete(key);
        else diversityState.hidden.add(key);
        renderInteractiveDiversityChart(chartWrap, seg, segData, genes, genomeLen);
        refreshLegend();
      });
    });
  };
  refreshLegend();
}

// What the quality filters did to this segment's genome set: how many survived, which were
// dropped and why, how many orientation-reversed genomes were fixed, which lineages were hidden for
// small n -- so a viewer can judge the curves instead of trusting them blindly.
function renderDiversityQc(segData) {
  const qc = segData.qc;
  if (!qc) return "";
  const list = (ids) => ids && ids.length
    ? `<details><summary>목록 보기 (${ids.length})</summary><div class="mono">${ids.map(escHtml).join(", ")}</div></details>` : "";
  const hidden = segData.lineages_hidden_small_n || {};
  const hiddenKeys = Object.keys(hidden);
  return `<div class="diversity-qc">
    <b>게놈 품질 필터</b> &mdash; 입력 ${qc.input_genomes}개 중 <b>${qc.kept}개</b> 사용
    (참조 대비 coverage &ge; ${qc.min_coverage}, identity &ge; ${qc.min_identity}).
    역방향(reverse)으로 들어온 게놈 <b>${qc.reversed_fixed}개</b>는 방향을 바로잡아 포함했습니다.
    <div>coverage 미달 ${(qc.excluded_low_coverage || []).length}개 ${list(qc.excluded_low_coverage)}</div>
    <div>identity 미달 ${(qc.excluded_low_identity || []).length}개 ${list(qc.excluded_low_identity)}</div>
    ${hiddenKeys.length ? `<div>n &lt; ${segData.min_lineage_n}라서 곡선을 숨긴 리니지: ${hiddenKeys.map((k) => `${escHtml(k)} (n=${hidden[k]})`).join(", ")} &mdash; 전체 LASV 곡선에는 포함됩니다.</div>` : ""}
  </div>`;
}

function renderDiversityPage() {
  const data = diversityState.data;
  if (!data) return "<p class='note'>불러오는 중...</p>";
  const availSegs = Object.keys(data.segments || {});
  if (!availSegs.length) {
    return `<h2>리니지별 유전체 다양성 (pixy &#960;)</h2>
      <p class="note">아직 diversity 데이터가 계산되지 않았습니다. Snakemake에서
      <code>output/diversity/{s,l}/pixy_pi.json</code>를 생성하는 <code>run_pixy</code> 규칙을
      실행한 후 다시 확인해주세요 (전체 파이프라인 실행 시 자동으로 포함됩니다).</p>`;
  }
  if (!diversityState.seg || !availSegs.includes(diversityState.seg)) diversityState.seg = availSegs[0];
  const seg = diversityState.seg;
  const segData = data.segments[seg];

  const segTabs = availSegs.map((s) => `<button class="seg-tab ${s === seg ? "active" : ""}" data-seg="${s}">${s.toUpperCase()} segment</button>`).join("");

  return `
    <h2>리니지별 유전체 다양성 (Nucleotide Diversity, pixy &#960;)</h2>
    <p class="note">Nextstrain의 genome-wide entropy 패널과 유사한 방식으로, population 게놈을 리니지별로
    묶어 <a href="https://github.com/ksamuk/pixy" target="_blank" rel="noopener">pixy</a>의 unbiased
    &#960; 추정치를 ${segData.window_size}bp 단위로 계산한 뒤 위 "평활화" 설정(기본 25bp 이동평균)으로 부드럽게 보여줍니다.
    참조 대비 coverage/identity가 낮은 게놈과 n이 작은 리니지는 제외했고(아래 품질 필터 참고), 역방향 게놈은 방향을 바로잡아 포함했습니다. 값이 높을수록 해당
    위치/리니지에서 유전적 다양성이 크다는 뜻이며, 프라이머/프로브 결합 부위가 다양성이 낮은 구간(보존
    영역)에 위치하는지 확인하는 데 참고할 수 있습니다.</p>
    <div class="seg-tab-row">${segTabs}</div>
    <div id="diversityChartWrap" class="diversity-chart-card"></div>
    <div id="diversityLegendWrap"></div>
    ${renderDiversityQc(segData)}
    <div id="diversityRegionAlignWrap"></div>
  `;
}

function wireDiversitySegTabs(view) {
  view.querySelectorAll(".seg-tab").forEach((btn) => {
    btn.addEventListener("click", () => {
      diversityState.seg = btn.dataset.seg;
      view.innerHTML = renderDiversityPage();
      wireDiversitySegTabs(view);
      mountDiversityChart(view);
    });
  });
}

async function showDiversityPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("diversityView");
  const view = $("#diversityView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    diversityState.data = await fetchJSON("/api/diversity");
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
    return;
  }
  view.innerHTML = renderDiversityPage();
  wireDiversitySegTabs(view);
  mountDiversityChart(view);
}

$("#diversityBtn").addEventListener("click", showDiversityPage);

// ---------------------------------------------------------------------------------------
// 유전자 보존영역: per gene of the reference genome, how conserved is it? Built from the same pixy pi tracks as the
// diversity page (one value per reference position, per lineage/clade and pooled over everything), so it needs no new
// pipeline step. A site is "변이 없음" when pi = 0 and "보존" when pi <= the threshold chosen on the page.
const gcSegTitle = (k) => (typeof segTitle === "function" ? segTitle(k) : `${String(k).toUpperCase()} segment`);
const gcLin = () => (typeof linKo === "function" ? linKo() : "리니지");   // the label differs per virus (Lineage / Clade)
const GC = { seg: null, thr: 0.01, minLen: 20, metric: "pi", sortKey: "start", sortDir: 1 };
const GC_GRADES = [
  { max: 0.01, label: "고도 보존", color: "#2e8b57" }, { max: 0.03, label: "보존", color: "#7fb069" },
  { max: 0.08, label: "중간", color: "#e0a030" }, { max: Infinity, label: "다양", color: "#d95f4b" },
];
const gcGrade = (pi) => GC_GRADES.find((g) => pi <= g.max);
const gcHeat = (v, metric) => {                // green = conserved, red = diverse
  if (v == null || Number.isNaN(v)) return "transparent";
  const t = metric === "pi" ? Math.min(1, v / 0.1) : 1 - Math.min(1, v);
  return `hsl(${Math.round(120 * (1 - t))}, 55%, 78%)`;
};
const gcFmt = (v, metric) => (v == null || Number.isNaN(v) ? "&mdash;" : metric === "pi" ? v.toFixed(4) : `${(v * 100).toFixed(1)}%`);

// position-indexed arrays {pi, w} (1-based) for every track of a segment
function gcTracks(segData, genomeLen) {
  const tracks = {};
  for (const [label, wins] of Object.entries(segData.windows || {})) {
    const pi = new Float32Array(genomeLen + 2).fill(NaN);
    for (const w of wins) {
      if (w.pi == null || !(w.no_sites > 0)) continue;
      for (let p = Math.max(1, w.start); p <= Math.min(genomeLen, w.end); p++) pi[p] = w.pi;
    }
    tracks[label] = pi;
  }
  return tracks;
}

// stats of one track over [start, end]: mean pi, share of invariant sites, share of conserved sites, longest conserved run
function gcStats(pi, start, end, thr) {
  let n = 0, sum = 0, inv = 0, cons = 0, run = 0, best = 0, bestEnd = 0;
  for (let p = start; p <= end; p++) {
    const v = pi[p];
    if (Number.isNaN(v)) { run = 0; continue; }
    n += 1; sum += v;
    if (v === 0) inv += 1;
    if (v <= thr) { cons += 1; run += 1; if (run > best) { best = run; bestEnd = p; } } else run = 0;
  }
  const len = end - start + 1;
  return { n, len, covered: len ? n / len : 0, mean: n ? sum / n : NaN, invariant: n ? inv / n : NaN, conserved: n ? cons / n : NaN,
    longest: best, longestStart: best ? bestEnd - best + 1 : 0, longestEnd: bestEnd };
}

// every conserved run (pi <= thr) of at least minLen positions
function gcRuns(pi, start, end, thr, minLen) {
  const runs = [];
  let s = null;
  for (let p = start; p <= end + 1; p++) {
    const ok = p <= end && !Number.isNaN(pi[p]) && pi[p] <= thr;
    if (ok && s === null) s = p;
    if (!ok && s !== null) { if (p - s >= minLen) runs.push({ start: s, end: p - 1, len: p - s }); s = null; }
  }
  return runs;
}

function gcRegions(data, seg) {
  const L = (data.genome_len || {})[seg] || 0;
  const genes = ((data.genes || {})[seg] || []).slice().sort((a, b) => a.start - b.start);
  const regions = genes.map((g) => ({ name: g.gene, start: g.start, end: g.end, gene: true }));
  const out = [];                                   // positions outside every gene (UTR / intergenic), as one row
  let cursor = 1;
  const gaps = [];
  for (const g of genes) { if (g.start > cursor) gaps.push([cursor, g.start - 1]); cursor = Math.max(cursor, g.end + 1); }
  if (cursor <= L) gaps.push([cursor, L]);
  return { L, regions, gaps, out };
}

function gcStatsAcross(tracks, ranges, thr) {      // ranges = [[start, end], ...] pooled into one region
  const res = {};
  for (const [label, pi] of Object.entries(tracks)) {
    let n = 0, sum = 0, inv = 0, cons = 0, len = 0;
    for (const [a, b] of ranges) {
      const s = gcStats(pi, a, b, thr);
      n += s.n; sum += s.mean * s.n || 0; inv += s.invariant * s.n || 0; cons += s.conserved * s.n || 0; len += s.len;
    }
    res[label] = { n, len, mean: n ? sum / n : NaN, invariant: n ? inv / n : NaN, conserved: n ? cons / n : NaN, longest: 0 };
  }
  return res;
}

function gcCompute() {
  const data = diversityState.data;
  const seg = GC.seg;
  const segData = data.segments[seg];
  const { L, regions, gaps } = gcRegions(data, seg);
  const tracks = gcTracks(segData, L);
  const pooledLabel = segData.all_label || PROFILE.population.pooled_label;
  const labels = Object.keys(tracks).filter((l) => l !== pooledLabel);
  const rows = regions.map((r) => ({ ...r, len: r.end - r.start + 1,
    stats: Object.fromEntries(Object.entries(tracks).map(([l, pi]) => [l, gcStats(pi, r.start, r.end, GC.thr)])) }));
  const outside = gaps.length ? { name: "유전자 밖 (UTR·유전자 사이)", start: null, end: null, gene: false,
    len: gaps.reduce((a, [s, e]) => a + (e - s + 1), 0), stats: gcStatsAcross(tracks, gaps, GC.thr) } : null;
  const whole = { name: "전체 유전체", start: 1, end: L, gene: false, len: L,
    stats: Object.fromEntries(Object.entries(tracks).map(([l, pi]) => [l, gcStats(pi, 1, L, GC.thr)])) };
  return { L, segData, tracks, pooledLabel, labels, rows, outside, whole };
}

async function gcTargetSets() {                    // which primer sets target which gene
  try {
    const o = await fetchJSON("/api/overview");
    const by = {};
    for (const s of o.sets || []) {
      const g = s.gene && s.gene.target_gene;
      if (g && g !== "NA") (by[g] = by[g] || []).push(s.display_name || s.set_id);
    }
    return by;
  } catch { return {}; }
}

function gcTrackSvg(c) {
  const W = 900, H = 150, padL = 36, padB = 40, plotW = W - padL - 10, plotH = H - padB - 8;
  const pooled = c.tracks[c.pooledLabel];
  const win = 25, smooth = new Float32Array(c.L + 2).fill(NaN);
  for (let p = 1; p <= c.L; p++) {
    let s = 0, n = 0;
    for (let q = Math.max(1, p - win); q <= Math.min(c.L, p + win); q++) if (!Number.isNaN(pooled[q])) { s += pooled[q]; n += 1; }
    smooth[p] = n ? s / n : NaN;
  }
  const ymax = Math.max(0.05, ...Array.from(smooth).filter((v) => !Number.isNaN(v)).map((v) => Math.min(v, 0.4)));
  const x = (p) => padL + ((p - 1) / Math.max(1, c.L - 1)) * plotW;
  const y = (v) => 8 + plotH * (1 - Math.min(v, ymax) / ymax);
  let d = "";
  for (let p = 1; p <= c.L; p += Math.max(1, Math.floor(c.L / 900))) if (!Number.isNaN(smooth[p])) d += `${d ? "L" : "M"}${x(p).toFixed(1)},${y(smooth[p]).toFixed(1)}`;
  const runs = gcRuns(pooled, 1, c.L, GC.thr, GC.minLen).map((r) => `<rect x="${x(r.start).toFixed(1)}" y="${8 + plotH + 2}" width="${Math.max(1.5, x(r.end) - x(r.start)).toFixed(1)}" height="6" fill="#2e8b57"><title>보존 구간 ${r.start}-${r.end} (${r.len} bp)</title></rect>`).join("");
  const genes = c.rows.map((r) => {
    const st = r.stats[c.pooledLabel], gr = gcGrade(st.mean);
    return `<g class="gc-gene" data-gene="${escHtml(r.name)}" style="cursor:pointer;"><rect x="${x(r.start).toFixed(1)}" y="${8 + plotH + 12}" width="${Math.max(2, x(r.end) - x(r.start)).toFixed(1)}" height="16" rx="3" fill="${gr.color}" opacity=".85">
      <title>${escHtml(r.name)} ${r.start}-${r.end} · 평균 π ${st.mean.toFixed(4)} · ${gr.label}</title></rect>
      <text x="${((x(r.start) + x(r.end)) / 2).toFixed(1)}" y="${8 + plotH + 24}" font-size="11" text-anchor="middle" fill="#fff" style="pointer-events:none;">${escHtml(r.name)}</text></g>`;
  }).join("");
  const ticks = [0, ymax / 2, ymax].map((v) => `<text x="${padL - 4}" y="${y(v) + 3}" font-size="9" text-anchor="end" fill="currentColor">${v.toFixed(2)}</text>
    <line x1="${padL}" x2="${W - 10}" y1="${y(v)}" y2="${y(v)}" stroke="currentColor" opacity=".12"></line>`).join("");
  return `<svg viewBox="0 0 ${W} ${H}" style="width:100%;color:var(--fg);">${ticks}<path d="${d}" fill="none" stroke="#4C78A8" stroke-width="1.4"></path>${runs}${genes}
    <text x="4" y="12" font-size="9" fill="currentColor">π</text></svg>
    <div class="note">파란 선: 전체(${escHtml(c.pooledLabel)}) π (±${win} 자리 평균) · 초록 막대: 연속 ${GC.minLen} bp 이상 보존 구간(π&le;${GC.thr}) · 아래 상자: 유전자(색 = 보존 등급, 눌러서 상세)</div>`;
}

function gcTable(c, targets) {
  const metric = GC.metric, pooled = c.pooledLabel;
  const val = (st) => (metric === "pi" ? st.mean : metric === "cons" ? st.conserved : st.invariant);
  const all = c.rows.slice();
  const keyOf = (r) => ({ start: r.start ?? 0, name: r.name, len: r.len, pooled: val(r.stats[pooled]) ?? 0,
    longest: r.stats[pooled].longest ?? 0 }[GC.sortKey] ?? (r.stats[GC.sortKey] ? val(r.stats[GC.sortKey]) : 0));
  all.sort((a, b) => { const x = keyOf(a), y = keyOf(b); return GC.sortDir * (typeof x === "string" ? x.localeCompare(y) : x - y); });
  const lineageHead = c.labels.map((l) => `<th data-gcsort="${escHtml(l)}" style="cursor:pointer;">${escHtml(l)}<div class="hm-total">n=${(c.segData.lineage_n || {})[l] ?? "?"}</div></th>`).join("");
  const cell = (st) => `<td style="background:${gcHeat(val(st), metric)};text-align:right;">${gcFmt(val(st), metric)}</td>`;
  const mkRow = (r, clickable) => {
    const p = r.stats[pooled], gr = gcGrade(p.mean);
    const tg = clickable && (targets[r.name] || []).length ? `<span class="pill pill-ok" title="${escHtml(targets[r.name].join(", "))}">프라이머 세트 ${targets[r.name].length}</span>` : "";
    return `<tr class="${clickable ? "gc-row" : ""}" data-gene="${escHtml(r.name)}" style="${clickable ? "cursor:pointer;" : ""}">
      <td><b>${escHtml(r.name)}</b> ${tg}</td><td class="mono">${r.start != null ? `${r.start}-${r.end}` : "&mdash;"}</td><td style="text-align:right;">${r.len.toLocaleString()}</td>
      ${cell(p)}<td style="text-align:right;">${p.longest ? `${p.longest} bp` : "&mdash;"}</td>
      <td><span class="pill" style="background:${gr.color};color:#fff;">${gr.label}</span></td>${c.labels.map((l) => cell(r.stats[l])).join("")}</tr>`;
  };
  return `<div style="overflow:auto;"><table class="mini-reference-table gc-table">
    <tr><th data-gcsort="name" style="cursor:pointer;">유전자</th><th data-gcsort="start" style="cursor:pointer;">위치</th><th data-gcsort="len" style="cursor:pointer;">길이</th>
      <th data-gcsort="pooled" style="cursor:pointer;">전체 ${escHtml(pooled)}<div class="hm-total">${metric === "pi" ? "평균 π" : metric === "cons" ? "보존 비율" : "변이 없는 자리"}</div></th>
      <th data-gcsort="longest" style="cursor:pointer;">최장 보존 구간</th><th>등급</th>${lineageHead}</tr>
    ${all.map((r) => mkRow(r, true)).join("")}${c.outside ? mkRow(c.outside, false) : ""}${mkRow(c.whole, false)}</table></div>`;
}

function gcRender(view, targets) {
  const data = diversityState.data;
  const segs = Object.keys(data.segments || {});
  if (!segs.length) { view.innerHTML = `<h2>유전자 보존영역</h2><p class="note">아직 다양성이 계산되지 않았습니다. "분석 다시하기"를 실행하세요.</p>`; return; }
  if (!GC.seg || !data.segments[GC.seg]) GC.seg = segs[0];
  const c = gcCompute();
  view._gc = { c, targets };
  view.innerHTML = `
    <h2>유전자 보존영역</h2>
    <p class="note">참조유전체의 유전자마다 ${gcLin()}별·전체에서 얼마나 보존되어 있는지 보여 줍니다. 값은 자리(염기)마다 계산한 nucleotide diversity(π)에서 나오며,
      π가 0이면 그 자리는 모든 게놈에서 같습니다. 낮을수록 보존, 높을수록 다양합니다. 보존된 연속 구간은 프라이머·프로브를 놓기 좋은 후보 자리입니다.</p>
    ${segs.length > 1 ? `<div class="seg-tab-row">${segs.map((s) => `<button class="seg-tab ${s === GC.seg ? "active" : ""}" data-gcseg="${s}">${escHtml(gcSegTitle(s))}</button>`).join("")}</div>` : ""}
    <div class="table-toolbar">
      <label>표시 지표 <select id="gcMetric"><option value="pi" ${GC.metric === "pi" ? "selected" : ""}>평균 π (낮을수록 보존)</option>
        <option value="cons" ${GC.metric === "cons" ? "selected" : ""}>보존 비율 (π &le; 기준인 자리)</option>
        <option value="inv" ${GC.metric === "inv" ? "selected" : ""}>변이 없는 자리 비율 (π = 0)</option></select></label>
      <label>보존 기준 π &le; <input id="gcThr" type="number" step="0.005" min="0" max="0.2" value="${GC.thr}" style="width:80px;"></label>
      <label>최소 보존 구간 <input id="gcMinLen" type="number" step="5" min="5" max="500" value="${GC.minLen}" style="width:70px;"> bp</label>
    </div>
    <div class="card">${gcTrackSvg(c)}</div>
    <div class="card"><h3>유전자별 보존 정도</h3>${gcTable(c, targets)}
      <div class="legend" style="margin-top:8px;">${GC_GRADES.map((g) => `<span><span class="dot" style="background:${g.color}"></span>${g.label} (전체 평균 π ${g.max === Infinity ? "&gt; 0.08" : "&le; " + g.max})</span>`).join("")}
        <span>열 제목을 누르면 정렬 · 행을 누르면 유전자 상세</span></div>
      <p class="note">${gcLin()} 열은 그 ${gcLin()} 안에서의 값입니다(게놈 수가 ${(c.segData.min_lineage_n ?? 10)}개 미만인 ${gcLin()}는 신뢰할 수 없어 제외). "전체"는 모든 게놈을 한 집단으로 본 값이라
        ${gcLin()} 사이의 차이까지 포함하므로 보통 개별 ${gcLin()}보다 높습니다.</p></div>`;
}

function gcWire(view) {
  view.addEventListener("click", (e) => {
    const seg = e.target.closest("[data-gcseg]"); if (seg) { GC.seg = seg.dataset.gcseg; gcRender(view, view._gc.targets); return; }
    const th = e.target.closest("[data-gcsort]");
    if (th) { const k = th.dataset.gcsort; GC.sortDir = GC.sortKey === k ? -GC.sortDir : (k === "name" || k === "start" ? 1 : 1); GC.sortKey = k; gcRender(view, view._gc.targets); return; }
    const g = e.target.closest("[data-gene]");
    if (g && g.dataset.gene) gcOpenGene(view._gc.c, view._gc.targets, g.dataset.gene);
  });
  view.addEventListener("change", (e) => {
    if (e.target.id === "gcMetric") GC.metric = e.target.value;
    else if (e.target.id === "gcThr") GC.thr = Math.max(0, Number(e.target.value) || 0);
    else if (e.target.id === "gcMinLen") GC.minLen = Math.max(5, Number(e.target.value) || 20);
    else return;
    gcRender(view, view._gc.targets);
  });
}

function gcOpenGene(c, targets, name) {
  const r = c.rows.find((x) => x.name === name);
  if (!r) return null;
  const palette = (typeof d3 !== "undefined" && d3.schemeTableau10) || ["#4C78A8", "#F2A65A", "#54A24B", "#E45756", "#B279A2", "#72B7B2"];
  const labels = [c.pooledLabel, ...c.labels];
  const rowsHtml = labels.map((l) => {
    const s = r.stats[l];
    return `<tr><td><b>${escHtml(l)}</b> <span class="note">n=${(c.segData.lineage_n || {})[l] ?? "?"}</span></td>
      <td style="background:${gcHeat(s.mean, "pi")};text-align:right;">${gcFmt(s.mean, "pi")}</td>
      <td style="background:${gcHeat(s.conserved, "cons")};text-align:right;">${gcFmt(s.conserved, "cons")}</td>
      <td style="background:${gcHeat(s.invariant, "inv")};text-align:right;">${gcFmt(s.invariant, "inv")}</td>
      <td style="text-align:right;">${s.longest ? `${s.longest} bp (${s.longestStart}-${s.longestEnd})` : "&mdash;"}</td></tr>`;
  }).join("");
  // per-lineage pi along the gene (smoothed)
  const W = 760, H = 170, pad = 34, len = r.end - r.start + 1, sm = Math.max(1, Math.round(len / 150));
  const series = labels.map((l) => {
    const pts = [];
    for (let p = r.start; p <= r.end; p += sm) {
      let s = 0, n = 0;
      for (let q = p; q < Math.min(r.end + 1, p + sm); q++) { const v = c.tracks[l][q]; if (!Number.isNaN(v)) { s += v; n += 1; } }
      pts.push([p, n ? s / n : NaN]);
    }
    return { l, pts };
  });
  const ymax = Math.max(0.05, ...series.flatMap((s) => s.pts.map((p) => p[1]).filter((v) => !Number.isNaN(v))));
  const x = (p) => pad + ((p - r.start) / Math.max(1, len - 1)) * (W - pad - 8), y = (v) => 8 + (H - 30) * (1 - Math.min(v, ymax) / ymax);
  const lines = series.map((s, i) => {
    const d = s.pts.filter((p) => !Number.isNaN(p[1])).map((p, k) => `${k ? "L" : "M"}${x(p[0]).toFixed(1)},${y(p[1]).toFixed(1)}`).join("");
    return `<path d="${d}" fill="none" stroke="${s.l === c.pooledLabel ? "currentColor" : palette[i % palette.length]}" stroke-width="${s.l === c.pooledLabel ? 2 : 1.3}" ${s.l === c.pooledLabel ? 'stroke-dasharray="5,3"' : ""}></path>`;
  }).join("");
  const legend = series.map((s, i) => `<span><span class="dot" style="background:${s.l === c.pooledLabel ? "var(--fg)" : palette[i % palette.length]}"></span>${escHtml(s.l)}</span>`).join("");
  const runs = gcRuns(c.tracks[c.pooledLabel], r.start, r.end, GC.thr, GC.minLen).sort((a, b) => b.len - a.len).slice(0, 8);
  const runRows = runs.map((u) => {
    const s = gcStats(c.tracks[c.pooledLabel], u.start, u.end, GC.thr);
    return `<tr><td class="mono">${u.start}-${u.end}</td><td style="text-align:right;">${u.len} bp</td><td style="text-align:right;">${s.mean.toFixed(4)}</td></tr>`;
  }).join("");
  const tg = targets[name] || [];
  return openDetailModal({
    title: `유전자 ${escHtml(name)} &mdash; 보존 상세`,
    subtitle: `${escHtml(c.segData.segment || GC.seg)} · ${r.start}-${r.end} (${len.toLocaleString()} bp) · 보존 기준 π&le;${GC.thr}`,
    html: `${tg.length ? `<p>이 유전자를 표적으로 하는 프라이머 세트: ${tg.map((n) => `<span class="pill pill-ok">${escHtml(n)}</span>`).join(" ")}</p>` : ""}
      <table class="mini-reference-table"><tr><th>${escHtml(gcLin())} / 집단</th><th>평균 π</th><th>보존 비율</th><th>변이 없는 자리</th><th>최장 보존 구간</th></tr>${rowsHtml}</table>
      <h4 style="margin:14px 0 4px;">유전자 위치별 π (${sm > 1 ? `${sm} bp 평균` : "자리별"})</h4>
      <svg viewBox="0 0 ${W} ${H}" style="width:100%;color:var(--fg);">${[0, ymax / 2, ymax].map((v) => `<line x1="${pad}" x2="${W - 8}" y1="${y(v)}" y2="${y(v)}" stroke="currentColor" opacity=".12"></line>
        <text x="${pad - 4}" y="${y(v) + 3}" font-size="9" text-anchor="end" fill="currentColor">${v.toFixed(2)}</text>`).join("")}${lines}
        <text x="${pad}" y="${H - 6}" font-size="10" fill="currentColor">${r.start}</text><text x="${W - 8}" y="${H - 6}" font-size="10" text-anchor="end" fill="currentColor">${r.end}</text></svg>
      <div class="legend">${legend}</div>
      <h4 style="margin:14px 0 4px;">가장 긴 보존 구간 (전체 기준, 최대 8개)</h4>
      ${runs.length ? `<table class="mini-reference-table"><tr><th>위치</th><th>길이</th><th>평균 π</th></tr>${runRows}</table>` : `<p class="note">이 유전자에는 연속 ${GC.minLen} bp 이상 π&le;${GC.thr}인 구간이 없습니다 (기준을 완화해 보세요).</p>`}`,
  });
}

async function showGeneConservationPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("geneConsView");
  const view = $("#geneConsView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    if (!diversityState.data) diversityState.data = await fetchJSON("/api/diversity");
    const targets = await gcTargetSets();
    if (!view._wired) { gcWire(view); view._wired = true; }
    gcRender(view, targets);
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`;
  }
}
$("#geneConsBtn").addEventListener("click", showGeneConservationPage);

// ---------------------------------------------------------------------------------------
// 리니지별 재설계 (varVAMP): per lineage (clade) a degenerate F/Probe/R scheme is designed from that lineage's alignment, then every
// design goes through the pipeline's real validation. Results exist only on this page (output/redesign/), not in the overview.
const RD = { seg: null, running: false, engine: "varvamp" };

function rdDesignRow(d, L, run) {
  const pct = (v) => (v == null ? "&mdash;" : `${Number(v).toFixed(1)}%`);
  const own = d.own_validated != null ? `<b>${d.own_validated} / ${d.own_total}</b> (${d.own_total ? (100 * d.own_validated / d.own_total).toFixed(1) : "0.0"}%)` : "&mdash;";
  const chips = (d.lineages || []).filter((x) => x.validated).map((x) => `<span class="pill ${x.lineage === L ? "pill-ok" : "pill-warn"}" title="${escHtml(x.lineage)}: ${x.validated} / ${x.total}">${escHtml(x.lineage)} ${x.validated}/${x.total}</span>`).join(" ");
  const seqCell = (k, label) => `<div><span class="note">${label}</span> <span class="mono">${escHtml(d[k])}</span> <button type="button" class="btn btn-secondary probe-copy" data-seq="${escHtml(d[k])}" style="padding:0 6px;" title="서열 복사">&#128203;</button> <span class="note" ${k === "P" && d.P_note ? `title="${escHtml(d.P_note)}"` : ""}>${k === "P" && d.P_note ? "⚠ " : ""}${d[k + "_tm"] != null ? `Tm ${d[k + "_tm"]}&deg;C` : ""}${d[k + "_gc"] != null ? ` · ${d[k + "_gc"]}%GC` : ""}${k === "P" && d.P_source ? ` · Probe: ${({ own: "자체 규칙", primer3: "Primer3", "own+primer3": "자체+Primer3" })[d.P_source] || d.P_source}` : ""}</span></div>`;
  return `<tr>
    <td class="mono">${escHtml(d.name)}</td>
    <td>${seqCell("F", "F")}${seqCell("P", "P")}${seqCell("R", "R")}</td>
    <td>${d.amplicon_len ?? "&mdash;"} bp<div class="note">${d.pair_cov != null ? `표본 쌍 커버리지 ${d.pair_cov}%` : `penalty ${d.penalty ?? "&mdash;"} · ΔG ${d.deltaG ?? "&mdash;"}`}</div></td>
    <td>${own}</td>
    <td>${d.validated ?? "&mdash;"} (${pct(d.validated_pct)})<div style="margin-top:2px;">${chips}</div></td>
    <td>${pct(d.probe_pct)}<div class="note">dropout ${d.dropout ?? "&mdash;"}</div></td>
    <td>${d.std_ok ? "✅" : "✗"}</td>
    <td><button type="button" class="btn btn-secondary rd-msa-btn" data-dir="${escHtml(d.dir)}" data-name="${escHtml(L + " / " + d.name)}">&#129516; MSA 보기</button></td></tr>`;
}

function renderRedesignResults(run, seg) {
  if (!run || !run.generated) return "<p class='note'>아직 실행하지 않았습니다. 위에서 리니지를 고르고 \"재설계 실행\"을 누르세요.</p>";
  const stale = state.hostScope && run.host_scope && state.hostScope !== run.host_scope;
  const head = `<p class="note">${escHtml(run.generated)} 실행 · <span class="pill pill-ok">${run.engine === "primerhunter" ? "PrimerHunter" : "varVAMP"}</span> ${pdScopeBadge(run.host_scope)} <span class="note">(집계 대상 ${run.population_n}개)</span>${run.done ? "" : " · <b>진행 중... 아래는 지금까지의 결과입니다</b>"}</p>`;
  const blocks = (run.lineages || []).map((e) => {
    const title = e.lineage === "ALL" ? "전체 (pan, 리니지 구분 없이)" : escHtml(e.lineage);
    if (e.status !== "ok") return `<div class="card" style="margin-top:10px;"><h4 style="margin:0;">${title} <span class="note">(${e.n}개)</span> <span class="pill pill-warn">${e.status === "skipped" ? "건너뜀" : "설계 없음"}</span></h4><p class="note">${escHtml(e.message || "")}</p></div>`;
    const cv = e.cover && e.cover.steps && e.cover.steps.length ? `<div style="margin:6px 0;padding:6px 10px;border:1px solid var(--border);border-radius:var(--r-md);background:rgba(84,162,74,.08);">
        <b>세트 커버</b> <span class="note">(설계를 함께 썼을 때 실제 검증을 통과하는 게놈 — 게놈을 가장 많이 추가하는 순서)</span>
        <div>${e.cover.steps.map((st, i) => `<span class="pill pill-ok" title="이 설계가 새로 추가한 게놈 ${st.added}개">${i ? "+ " : ""}${escHtml(st.name)} <b>+${st.added}</b> → ${st.covered} (${st.pct}%)</span>`).join(" ")}</div>
        <div style="margin-top:4px;"><b>합계 ${e.cover.covered} / ${e.cover.total} (${e.cover.pct}%)</b> <span class="note">설계 ${e.cover.steps.length}개 사용 · 위 표의 "이 리니지 Validated"는 설계 하나씩의 값입니다</span></div></div>` : "";
    return `<div class="card" style="margin-top:10px;"><h4 style="margin:0;">${title} <span class="note">(${e.n}개 게놈 · ${e.sample_n ? `PrimerHunter 대표 표본 ${e.sample_n}개` : `consensus 임계값 ${e.threshold}`})</span></h4>${cv}
      <div style="overflow:auto;"><table class="mini-reference-table"><tr><th>설계</th><th>서열 (5′→3′, IUPAC 포함)</th><th>Amplicon</th><th>${e.lineage === "ALL" ? "" : "이 리니지 Validated"}</th><th>전체 Validated <span class="note">(통과한 리니지)</span></th><th>Probe 일치</th><th>표준물질</th><th></th></tr>
      ${e.designs.map((d) => rdDesignRow(d, e.lineage, run)).join("")}</table></div></div>`;
  }).join("");
  return head + (stale ? `<p class="note"><span class="pill pill-warn">지금 화면의 숙주 기준과 다릅니다</span> 이 결과는 ${escHtml(PD_SCOPE_LABEL[run.host_scope] || run.host_scope)} 기준입니다.</p>` : "") + blocks;
}

async function showRedesignPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("redesignView");
  const view = $("#redesignView");
  const segs = Object.keys((state.hostInfo || {}).segments || {}).map((x) => x.toUpperCase());
  const allSegs = segs.length ? segs : ["S", "L"];
  if (!RD.seg || !allSegs.includes(RD.seg)) RD.seg = allSegs[0];
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  let run;
  try { run = await fetchJSON(`/api/redesign/${RD.seg.toLowerCase()}/run?engine=${RD.engine}`); } catch (e) { view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`; return; }
  const avail = run.available || {};
  const checked = RD.checked || {};
  view.innerHTML = `
    <h2>🛠 리니지별 재설계</h2>
    <p class="note">리니지(클레이드)마다 <b>varVAMP</b>(그 리니지의 전체 게놈 정렬, qPCR 모드로 degenerate F/R + Probe) 또는 <b>PrimerHunter</b>(대표 게놈 표본에서 Tm 기반으로 F/R 쌍과 쌍 집합을 설계하고 Probe는 이 도구의 Probe 설계기로 추가)로 프라이머·Probe를 새로 설계하고, 설계된 세트를 <b>이 파이프라인의 실제 검증</b>(F/R mismatch ≤ 기준 + 3′ 끝, Probe 완전 일치, 표준물질, dimer/hairpin)에 태웁니다.
      consensus 임계값은 높은 값부터 낮춰 가며 설계가 나올 때까지 시도합니다. 결과는 이 페이지에서만 보이며 개요·요약·보고서에는 들어가지 않습니다. 설계는 후보일 뿐이라 합성 전 실험 확인이 필요합니다.</p>
    <div class="seg-tab-row">${allSegs.map((x) => `<button class="seg-tab ${x === RD.seg ? "active" : ""}" data-rdseg="${x}">${x}</button>`).join("")}</div>
    <div class="card">
      <div class="table-toolbar" style="flex-wrap:wrap;gap:12px;">
        <span class="note">숙주 기준 ${pdScopeBadge(run.active_scope)}</span>
        <label class="note" title="설계 도구. PrimerHunter는 느립니다(리니지당 수 분~수십 분)">도구 <select id="rdEngine"><option value="varvamp" ${RD.engine === "varvamp" ? "selected" : ""}>varVAMP</option><option value="primerhunter" ${RD.engine === "primerhunter" ? "selected" : ""}>PrimerHunter</option></select></label>
        <label class="note" ${RD.engine === "primerhunter" ? "" : 'hidden'} title="PrimerHunter는 모든 게놈이 아니라 가장 다양한 대표 게놈 N개로 설계합니다 (검증은 리니지 전체로)">대표 표본 <input id="rdSample" type="number" min="5" max="100" value="${RD.sample ?? 25}" style="width:56px;"></label>
        <label class="note" ${RD.engine === "primerhunter" ? "" : 'hidden'} title="프라이머-표본 결합의 최소 Tm(℃). 높일수록 mismatch가 적은(우리 검증 기준에 가까운) 프라이머만 남지만 후보가 줄어듭니다">최소 Tm <input id="rdMinTm" type="number" min="40" max="70" value="${RD.minTm ?? 56}" style="width:56px;"></label>
        <label class="note" ${RD.engine === "primerhunter" ? "" : 'hidden'} title="한 프라이머가 이 비율 이상의 표본에 결합해야 합니다. 100 미만이면 여러 쌍으로 전체를 덮는 세트 커버를 씁니다">최소 커버리지 % <input id="rdMinCov" type="number" min="10" max="100" value="${RD.minCov ?? 80}" style="width:56px;"></label>
        <label class="note" title="consensus 임계값: 높을수록 특이적이지만 설계가 안 나올 수 있습니다">임계값 <input id="rdThr" type="number" step="0.05" min="0.3" max="1" value="${RD.thr ?? 0.9}" style="width:62px;"> → 최저 <input id="rdMinThr" type="number" step="0.05" min="0.3" max="1" value="${RD.minThr ?? 0.6}" style="width:62px;"></label>
        <label class="note" title="F/R 프라이머 하나에 넣을 수 있는 모호 염기(IUPAC) 최대 수">F/R IUPAC <input id="rdAmbig" type="number" min="0" max="8" value="${RD.ambig ?? 2}" style="width:48px;"></label>
        <label class="note" title="Probe의 모호 염기 최대 수 (프라이머보다 엄격하게 1 적게가 권장)">Probe IUPAC <input id="rdPAmbig" type="number" min="0" max="6" value="${RD.pAmbig ?? 1}" style="width:48px;"></label>
        <label class="note" title="리니지마다 검증할 설계 수">설계 수 <input id="rdTop" type="number" min="1" max="8" value="${RD.top ?? 3}" style="width:48px;"></label>
        <label class="note" title="게놈이 이보다 적은 리니지는 건너뜁니다">최소 게놈 수 <input id="rdMinN" type="number" min="3" max="1000" value="${RD.minN ?? 5}" style="width:56px;"></label>
        <button id="rdRunBtn" class="btn btn-primary">&#9654; 재설계 실행</button>
      </div>
      <div style="margin-top:8px;"><b class="note">리니지 선택</b>
        ${Object.entries(avail).map(([L, n]) => `<label class="note" style="margin-right:12px;"><input type="checkbox" class="rd-lin" value="${escHtml(L)}" ${(checked[L] ?? n >= (RD.minN ?? 5)) ? "checked" : ""}> ${escHtml(L)} <span class="mono">${n}</span></label>`).join("")}
        <label class="note" style="margin-right:12px;"><input type="checkbox" id="rdAll" ${RD.all ? "checked" : ""}> 전체(pan)도 함께</label></div>
      <div id="rdStatus" class="note" style="margin-top:6px;"></div>
    </div>
    <div id="rdResults">${renderRedesignResults(run, RD.seg)}</div>`;
  const results = view.querySelector("#rdResults"), statusEl = view.querySelector("#rdStatus"), runBtn = view.querySelector("#rdRunBtn");
  view.querySelectorAll("[data-rdseg]").forEach((b) => b.addEventListener("click", () => { RD.seg = b.dataset.rdseg; showRedesignPage(); }));
  view.querySelector("#rdEngine").addEventListener("change", (e) => { RD.engine = e.target.value; showRedesignPage(); });
  results.addEventListener("click", (e) => {
    const c = e.target.closest(".probe-copy");
    if (c) { copyTextToClipboard(c.dataset.seq).then(() => { c.textContent = "✓"; setTimeout(() => { c.innerHTML = "&#128203;"; }, 1000); }).catch(() => uiAlert("복사하지 못했습니다.")); return; }
    const m = e.target.closest(".rd-msa-btn");
    if (m) openProbeMsaModal({ set_id: m.dataset.name }, "", `/api/redesign/${RD.seg.toLowerCase()}/msa?engine=${RD.engine}&host=${encodeURIComponent(state.hostScope || "")}&dir=${encodeURIComponent(m.dataset.dir)}`,
      `전체 게놈 정렬 (MSA) &mdash; 재설계 ${escHtml(m.dataset.name)}`);
  });
  const poll = async (jobId) => {
    RD.running = true; runBtn.disabled = true;
    const jt = jobToast("재설계");
    for (;;) {
      let j;
      try { j = await fetchJSON(`/api/jobs/${jobId}`); } catch (e) { break; }
      if (!document.body.contains(results)) { RD.running = false; jt.done(j.status !== "error"); return; }
      if (j.status !== "running") { if (j.status === "error") uiAlert("재설계에 실패했습니다: " + (j.error || "")); jt.done(j.status === "done"); break; }
      jt.progress(j);
      statusEl.textContent = `⏳ ${j.progress || ""} (${j.elapsed_sec ?? 0}초)`;
      try { results.innerHTML = renderRedesignResults(await fetchJSON(`/api/redesign/${RD.seg.toLowerCase()}/run?engine=${RD.engine}`), RD.seg); } catch (e) { /* keep the last view */ }
      await new Promise((r) => setTimeout(r, 3000));
    }
    RD.running = false; runBtn.disabled = false; statusEl.textContent = "";
    if (document.body.contains(results)) results.innerHTML = renderRedesignResults(await fetchJSON(`/api/redesign/${RD.seg.toLowerCase()}/run?engine=${RD.engine}`), RD.seg);
  };
  runBtn.addEventListener("click", async () => {
    const num = (id) => Number(view.querySelector(id).value);
    Object.assign(RD, { thr: num("#rdThr"), minThr: num("#rdMinThr"), ambig: num("#rdAmbig"), pAmbig: num("#rdPAmbig"), top: num("#rdTop"), minN: num("#rdMinN"), all: view.querySelector("#rdAll").checked, sample: num("#rdSample"), minCov: num("#rdMinCov"), minTm: num("#rdMinTm") });
    RD.checked = Object.fromEntries([...view.querySelectorAll(".rd-lin")].map((x) => [x.value, x.checked]));
    const lineages = [...view.querySelectorAll(".rd-lin:checked")].map((x) => x.value);
    if (!lineages.length && !RD.all) { uiAlert("리니지를 하나 이상 고르세요."); return; }
    if (!(await uiConfirm(`${lineages.length}개 리니지${RD.all ? " + 전체" : ""}를 설계하고 설계마다 실제 검증을 돌립니다. ${RD.engine === "primerhunter" ? "PrimerHunter는 리니지당 수 분에서 수십 분이 걸리고 CPU를 여러 개 씁니다." : "리니지당 1분 안팎이 걸립니다."} 그동안 서버가 계산에 쓰입니다. 계속할까요?`, { title: "리니지별 재설계", okLabel: "시작" }))) return;
    try {
      const { job_id, error } = await fetchJSON("/api/redesign", { method: "POST", headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ segment: RD.seg.toLowerCase(), host: state.hostScope, lineages, include_all: RD.all, threshold: RD.thr, min_threshold: RD.minThr, ambig: RD.ambig, probe_ambig: RD.pAmbig, top: RD.top, min_n: RD.minN, engine: RD.engine, ph_sample: RD.sample, ph_min_cov: RD.minCov, ph_min_tm: RD.minTm }) });
      if (error) throw new Error(error);
      poll(job_id);
    } catch (err) { uiAlert("재설계를 시작하지 못했습니다: " + err.message); }
  });
  try {                                                         // pick up a run that is already going (started from another computer)
    const active = (await fetchJSON("/api/jobs/active")).find((j) => j.kind === "redesign");
    if (active) poll(active.job_id);
  } catch (e) { /* ignore */ }
}
$("#redesignBtn").addEventListener("click", showRedesignPage);

// 사이드바의 메뉴 버튼 전체 접기/펼치기(프라이머 그룹·세트 목록은 그대로 보임): 이 브라우저에 기억합니다.
(function initMenuToggle() {
  const btn = document.getElementById("menuToggleBtn"), menu = document.getElementById("sidebarMenu");
  if (!btn || !menu) return;
  const apply = (hidden) => {
    menu.hidden = hidden;
    btn.innerHTML = hidden ? "&#9776; 메뉴 펼치기 &#9662;" : "&#9776; 메뉴 접기 &#9652;";
    try { localStorage.setItem("lasv_menu_hidden", hidden ? "1" : "0"); } catch (e) { /* ignore */ }
  };
  btn.addEventListener("click", () => apply(!menu.hidden));
  let saved = null;
  try { saved = localStorage.getItem("lasv_menu_hidden"); } catch (e) { /* ignore */ }
  apply(saved === "1");
})();

// 사이드바 메뉴 묶음(분석·설계·정보·관리): 열림/닫힘을 이 브라우저에 기억합니다.
(function initMenuGroups() {
  let saved = {};
  try { saved = JSON.parse(localStorage.getItem("lasv_menu_open") || "{}"); } catch (e) { /* ignore */ }
  document.querySelectorAll(".menu-group").forEach((d) => {
    const k = d.dataset.menu;
    if (k in saved) d.open = !!saved[k];
    d.addEventListener("toggle", () => {
      saved[k] = d.open;
      try { localStorage.setItem("lasv_menu_open", JSON.stringify(saved)); } catch (e) { /* ignore */ }
    });
  });
})();

// ---------------------------------------------------------------------------------------
// 사이드바 접기 / 펼치기 (« / », Ctrl+B) -- 선택은 이 브라우저에 기억됩니다.
(function initSidebarCollapse() {
  const app = $("#appRoot"), openBtn = $("#sidebarOpenBtn"), closeBtn = $("#sidebarCollapseBtn"), edgeBtn = $("#sidebarEdgeBtn");
  let expandedW = 320;
  const placeEdge = (collapsed) => {                  // the handle sits on the sidebar's right edge, whatever width the sidebar has
    const sb = document.querySelector(".sidebar");
    if (!collapsed && sb.offsetWidth > 40) expandedW = sb.offsetWidth;
    edgeBtn.style.left = collapsed ? "0px" : `${expandedW - 1}px`;
  };
  window.addEventListener("resize", () => placeEdge(app.classList.contains("sidebar-collapsed")));
  const apply = (collapsed) => {
    app.classList.toggle("sidebar-collapsed", collapsed);
    openBtn.hidden = !collapsed;
    edgeBtn.innerHTML = collapsed ? "&rsaquo;" : "&lsaquo;";
    placeEdge(collapsed);
    try { localStorage.setItem("lasv_sidebar_collapsed", collapsed ? "1" : "0"); } catch (e) { /* private mode */ }
  };
  closeBtn.addEventListener("click", () => apply(true));
  edgeBtn.addEventListener("click", () => apply(!app.classList.contains("sidebar-collapsed")));
  openBtn.addEventListener("click", () => apply(false));
  document.addEventListener("keydown", (e) => {
    if ((e.ctrlKey || e.metaKey) && !e.shiftKey && !e.altKey && e.key.toLowerCase() === "b" && !/^(INPUT|TEXTAREA|SELECT)$/.test((document.activeElement || {}).tagName || "")) {
      e.preventDefault(); apply(!app.classList.contains("sidebar-collapsed"));
    }
  });
  let saved = null;
  try { saved = localStorage.getItem("lasv_sidebar_collapsed"); } catch (e) { /* ignore */ }
  apply(saved === "1");
})();

// parent / variants / identical-sequence facts of a set (computed by the server, see annotate_set_relations)
function relationsHtml(setId) {
  const me = state.sets.find((x) => x.set_id === setId);
  if (!me) return "";
  const nameOf = (id) => ((state.sets.find((x) => x.set_id === id) || {}).display_name || id);
  const link = (id) => `<a href="#" class="rel-link" data-set-id="${escHtml(id)}">${escHtml(nameOf(id))}</a>`;
  const parts = [];
  if (me.parent_id) parts.push(`<b>변이 세트</b> · 부모: ${link(me.parent_id)} <span class="note">(${me.parent_origin === "manual" ? "직접 지정" : "이름 규칙 \"-숫자\"로 추정"})</span>`);
  const variants = state.sets.filter((x) => x.parent_id === setId);
  if (variants.length) parts.push(`<b>변이 ${variants.length}개</b>: ${variants.map((v) => link(v.set_id)).join(" · ")}`);
  if (me.duplicate_ids && me.duplicate_ids.length) parts.push(`<span class="dup-badge">동일 서열</span> F/R/Probe가 같은 세트: ${me.duplicate_ids.map(link).join(" · ")} <span class="note">(이름만 다른 중복 등록이면 하나만 남기거나 그룹에서 한쪽만 쓰세요)</span>`);
  return `<div class="relations-box note">${parts.map((x) => `<div>${x}</div>`).join("")}<div><button type="button" class="btn btn-secondary rel-parent-edit" data-set-id="${escHtml(setId)}" style="padding:1px 10px;">부모 세트 ${me.parent_id ? "변경" : "지정"}</button></div></div>`;
}

document.addEventListener("click", async (e) => {
  const lk = e.target.closest(".rel-link");
  if (lk) { e.preventDefault(); selectSet(lk.dataset.setId); return; }
  const eb = e.target.closest(".rel-parent-edit");
  if (!eb) return;
  const me = state.sets.find((x) => x.set_id === eb.dataset.setId);
  if (!me) return;
  const cands = state.sets.filter((x) => x.segment === me.segment && x.set_id !== me.set_id && !x.parent_id)
    .sort((a, b) => (a.display_name || a.set_id).localeCompare(b.display_name || b.set_id, undefined, { numeric: true }));
  const { close, overlay } = openDetailModal({
    title: "부모 세트 지정", subtitle: `${escHtml(me.display_name || me.set_id)}가 어느 세트에서 변형한 것인지 고릅니다 (한 단계만 지원)`, compact: true,
    html: `<select id="relParentSel" style="width:100%;margin:8px 0;"><option value="">(부모 없음)</option>${cands.map((c) => `<option value="${escHtml(c.set_id)}" ${c.set_id === me.parent_id ? "selected" : ""}>${escHtml(c.display_name || c.set_id)}</option>`).join("")}</select>
      <button type="button" class="btn btn-primary" id="relParentSave">저장</button> <span class="note" id="relParentMsg"></span>`,
  });
  overlay.querySelector("#relParentSave").addEventListener("click", async () => {
    const parent = overlay.querySelector("#relParentSel").value;
    try {
      const r = await fetch(`/api/set-parent/${encodeURIComponent(me.set_id)}`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ parent }) });
      const out = await r.json();
      if (!r.ok) throw new Error(out.error || r.status);
      close(); await loadSets(); selectSet(me.set_id);
    } catch (err) { overlay.querySelector("#relParentMsg").textContent = "저장 실패: " + err.message; }
  });
});

// ---------------------------------------------------------------------------------------
// 사용자의 최종 결과(판정): 세트·그룹마다 채택 / 보류 / 탈락 + 한 줄 메모. 이 대시보드를 쓰는 사람이 직접 정하는 표시이며 파이프라인 계산에는 쓰이지 않습니다.
const FR_LABEL = { adopt: "✅ 채택", hold: "⏸ 보류", reject: "✖ 탈락" };
state.finalResults = {};
async function loadFinalResults() {
  try { state.finalResults = await fetchJSON("/api/final-results"); } catch (e) { /* keep the last known */ }
}
function finalChipHtml(id) {
  const r = (state.finalResults || {})[id];
  if (!r || (!r.verdict && !r.note)) return '<span class="note">&mdash;</span>';
  return `${r.verdict ? `<span class="fr-chip ${r.verdict}" title="${escHtml(r.updated || "")}">${FR_LABEL[r.verdict]}</span>` : ""}${r.note ? `<div class="note" title="${escHtml(r.note)}" style="max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;">${escHtml(r.note)}</div>` : ""}`;
}
function finalResultBoxHtml(id) { return `<div class="final-result-box" data-fr-id="${escHtml(id)}"></div>`; }
function initFinalResult(root) {
  const box = root && root.querySelector(".final-result-box");
  if (!box) return;
  const id = box.dataset.frId;
  const draw = () => {
    const r = (state.finalResults || {})[id] || { verdict: "", note: "" };
    box.innerHTML = `<span class="fr-label">최종 결과</span>
      ${["adopt", "hold", "reject"].map((v) => `<button type="button" class="fr-btn ${v} ${r.verdict === v ? "on" : ""}" data-verdict="${v}">${FR_LABEL[v]}</button>`).join("")}
      <button type="button" class="fr-btn" data-verdict="" title="판정을 지웁니다">미정</button>
      <input type="text" class="fr-note" maxlength="500" placeholder="사유·메모 (예: 실험 후 채택, GPC 쪽이 더 안정적...)" value="${escHtml(r.note || "")}">
      <span class="note fr-status">${r.updated ? `마지막 저장 ${escHtml(r.updated)}` : "아직 표시하지 않았습니다"}</span>`;
  };
  const save = async (verdict, note) => {
    const st = box.querySelector(".fr-status");
    if (st) st.textContent = "저장 중...";
    try {
      const out = await fetchJSON(`/api/final-results/${encodeURIComponent(id)}`, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ verdict, note }) });
      if (out.result) state.finalResults[id] = out.result; else delete state.finalResults[id];
      draw();
    } catch (e) { if (st) st.textContent = "저장 실패: " + e.message; }
  };
  draw();
  box.addEventListener("click", (e) => {
    const b = e.target.closest(".fr-btn");
    if (b) save(b.dataset.verdict, box.querySelector(".fr-note").value);
  });
  box.addEventListener("change", (e) => {
    if (e.target.classList.contains("fr-note")) save(((state.finalResults || {})[id] || {}).verdict || "", e.target.value);
  });
}

// ---------------------------------------------------------------------------------------
// 커버리지 벤 다이어그램: which genomes do the selected primer sets / groups validate -- union, intersection, every overlap region and the genomes
// no selected set catches (population minus union). Built from /api/venn (validated ids per set and group, in the active host scope).
const VN = { seg: null, selected: null, lineage: "ALL", data: null, region: null };
const VN_COLORS = ["#4C78A8", "#F2A65A", "#54A24B", "#B279A2", "#E45756", "#72B7B2", "#9D755D", "#BAB0AC"];

function vnEllipses(k) {
  // unit-space geometry [cx, cy, a, b, rotation deg]; scale/offset to pixels in vnToPx
  if (k === 1) return { es: [[0, 0, 1, 1, 0]], scale: 130 };
  if (k === 2) return { es: [[-0.55, 0, 1, 1, 0], [0.55, 0, 1, 1, 0]], scale: 125 };
  if (k === 3) return { es: [[-0.55, -0.35, 1, 1, 0], [0.55, -0.35, 1, 1, 0], [0, 0.6, 1, 1, 0]], scale: 105 };
  return { es: [[-0.7, -0.3, 2.0, 1.1, 45], [-0.2, -0.7, 2.0, 1.1, 45], [0.2, -0.7, 2.0, 1.1, -45], [0.7, -0.3, 2.0, 1.1, -45]], scale: 78 };
}
function vnInside(e, x, y) {
  const t = (e[4] * Math.PI) / 180, dx = x - e[0], dy = y - e[1];
  const u = dx * Math.cos(t) + dy * Math.sin(t), w = -dx * Math.sin(t) + dy * Math.cos(t);
  return (u / e[2]) ** 2 + (w / e[3]) ** 2 <= 1;
}
function vnMask(es, x, y) { let m = 0; es.forEach((e, i) => { if (vnInside(e, x, y)) m |= 1 << i; }); return m; }

function vnCompute() {
  const d = VN.data, byKey = new Map(d.items.map((it) => [it.key, it]));
  const sel = [...VN.selected].map((k) => byKey.get(k)).filter(Boolean);          // selection order = colour order, everywhere on the page
  const inLin = (lin) => VN.lineage === "ALL" || lin === VN.lineage;
  const sets = sel.map((it) => new Set(it.ids));
  const maskOf = (id) => { let m = 0; sets.forEach((st, i) => { if (st.has(id)) m |= 1 << i; }); return m; };
  const byMask = new Map(), all = [];
  const perLin = new Map();
  for (const [id, lin, country] of d.genomes) {
    const m = maskOf(id);
    const rec = { id, lin, country, m };
    if (!perLin.has(lin)) perLin.set(lin, { n: 0, union: 0, inter: 0, none: 0, each: sel.map(() => 0) });
    const pl = perLin.get(lin);
    pl.n += 1; if (m) pl.union += 1; else pl.none += 1; if (sel.length && m === (1 << sel.length) - 1) pl.inter += 1;
    sel.forEach((_, i) => { if (m & (1 << i)) pl.each[i] += 1; });
    if (!inLin(lin)) continue;
    all.push(rec);
    if (!byMask.has(m)) byMask.set(m, []);
    byMask.get(m).push(rec);
  }
  const n = all.length, full = (1 << sel.length) - 1;
  const union = all.filter((r) => r.m).length, inter = sel.length ? (byMask.get(full) || []).length : 0;
  return { sel, byMask, all, n, union, inter, none: n - union, full, perLin };
}

// the circles' colours, named: one chip per selected set/group in the same colour as its circle, its column in the lineage table and its checkbox
function vnLegendHtml(c) {
  return `<div class="vn-legend">${c.sel.map((x, i) => `<span class="vn-chip"><i style="background:${VN_COLORS[i % VN_COLORS.length]}"></i>${escHtml(x.name)} <b>${x.ids.length}</b></span>`).join("")}</div>`;
}

function vnRegionName(m, sel) {
  if (!m) return "어느 세트도 못 잡음";
  const names = sel.filter((_, i) => m & (1 << i)).map((x) => x.name);
  const only = names.length === 1;
  return only ? `${names[0]} 단독` : names.join(" ∩ ");
}

function vnDiagram(c) {
  const k = c.sel.length;
  if (k < 1 || k > 4) return "";
  const { es, scale } = vnEllipses(k), W = 560, H = 430, cx = 280, cy = 215;
  const toPx = (x, y) => [cx + x * scale, cy + y * scale];
  // a representative point inside every region (the sampled point of the region nearest to its centroid)
  const pts = new Map();
  for (let i = -330; i <= 330; i += 3) for (let j = -330; j <= 330; j += 3) {
    const x = i / 100, y = j / 100, m = vnMask(es, x, y);
    if (!m) continue;
    if (!pts.has(m)) pts.set(m, []);
    pts.get(m).push([x, y]);
  }
  const anchor = (m) => {
    const a = pts.get(m); if (!a || !a.length) return null;
    const mx = a.reduce((s2, p) => s2 + p[0], 0) / a.length, my = a.reduce((s2, p) => s2 + p[1], 0) / a.length;
    return a.reduce((best, p) => ((p[0] - mx) ** 2 + (p[1] - my) ** 2 < (best[0] - mx) ** 2 + (best[1] - my) ** 2 ? p : best), a[0]);
  };
  const shapes = es.map((e, i) => {
    const [px, py] = toPx(e[0], e[1]);
    return `<ellipse cx="${px}" cy="${py}" rx="${e[2] * scale}" ry="${e[3] * scale}" transform="rotate(${e[4]} ${px} ${py})" fill="${VN_COLORS[i]}" fill-opacity=".22" stroke="${VN_COLORS[i]}" stroke-width="2.5"/>`;
  }).join("");
  const labels = [...pts.keys()].map((m) => {
    const a = anchor(m), cnt = (c.byMask.get(m) || []).length, [px, py] = toPx(a[0], a[1]);
    const pct = c.n ? (100 * cnt / c.n).toFixed(1) : "0.0";
    const hot = VN.region === m;
    return `<g class="vn-region-label" data-mask="${m}" style="cursor:pointer;"><title>${escHtml(vnRegionName(m, c.sel))}: ${cnt}개 (${pct}%)</title>
      <text x="${px}" y="${py}" text-anchor="middle" font-size="${k === 4 ? 14 : 18}" font-weight="700" fill="currentColor"${hot ? ' text-decoration="underline"' : ""}>${cnt}</text>
      ${k <= 3 ? `<text x="${px}" y="${py + 15}" text-anchor="middle" font-size="11" fill="currentColor" opacity=".7">${pct}%</text>` : ""}</g>`;
  }).join("");
  const none = c.none, nonePct = c.n ? (100 * none / c.n).toFixed(1) : "0.0";
  return `<svg id="vnSvg" viewBox="0 0 ${W} ${H}" style="width:100%;max-width:640px;height:auto;display:block;margin:auto;color:var(--fg);">
    <rect x="4" y="4" width="${W - 8}" height="${H - 8}" rx="12" fill="none" stroke="var(--border)" stroke-width="2"/>
    <text x="16" y="24" font-size="12" fill="currentColor" opacity=".75">전체 ${c.n}개${VN.lineage === "ALL" ? "" : ` (${escHtml(VN.lineage)})`}</text>
    ${shapes}${labels}
    <g class="vn-region-label" data-mask="0" style="cursor:pointer;"><title>어느 세트도 못 잡는 유전체</title>
      <rect x="${W - 150}" y="${H - 52}" width="138" height="40" rx="8" fill="var(--bad)" fill-opacity=".14" stroke="var(--bad)"/>
      <text x="${W - 81}" y="${H - 34}" text-anchor="middle" font-size="12" fill="currentColor">못 잡음 (여집합)</text>
      <text x="${W - 81}" y="${H - 18}" text-anchor="middle" font-size="14" font-weight="700" fill="currentColor">${none} (${nonePct}%)</text></g>
  </svg>`;
}

function vnRegionTable(c) {
  const rows = [...c.byMask.entries()].sort((a, b) => b[1].length - a[1].length).map(([m, arr]) => {
    const pct = c.n ? (100 * arr.length / c.n).toFixed(1) : "0.0";
    const bad = m === 0;
    return `<tr class="vn-region-row" data-mask="${m}" style="cursor:pointer;${bad ? "background:rgba(228,87,86,.08);" : ""}"><td>${escHtml(vnRegionName(m, c.sel))}</td>
      <td>${c.sel.map((_, i) => `<span class="dot" style="background:${m & (1 << i) ? VN_COLORS[i] : "transparent"};border:1px solid ${VN_COLORS[i]};display:inline-block;width:10px;height:10px;border-radius:50%;margin-right:2px;"></span>`).join("")}</td>
      <td><b>${arr.length}</b></td><td>${pct}%</td></tr>`;
  }).join("");
  return `<table class="mini-reference-table"><tr><th>영역 (겹침 조합)</th><th></th><th>유전체 수</th><th>전체 대비</th></tr>${rows}</table>`;
}

function vnLineageTable(c) {
  const lins = [...c.perLin.keys()].sort((a, b) => a.localeCompare(b));
  const pct = (v, n) => (n ? (100 * v / n).toFixed(1) : "0.0");
  const head = `<tr><th>${escHtml((PROFILE.lineage && PROFILE.lineage.label) || "Lineage")}</th><th>유전체</th>${c.sel.map((x, i) => `<th style="border-bottom:3px solid ${VN_COLORS[i]};">${escHtml(x.name)}</th>`).join("")}<th>합집합</th><th>교집합</th><th>못 잡음</th></tr>`;
  const body = lins.map((L) => {
    const r = c.perLin.get(L), cur = VN.lineage === L;
    return `<tr${cur ? ' style="background:rgba(76,120,168,.12);"' : ""}><td><a href="#" class="vn-pick-lin" data-lin="${escHtml(L)}">${escHtml(L)}</a></td><td>${r.n}</td>
      ${r.each.map((v) => `<td>${v} <span class="note">${pct(v, r.n)}%</span></td>`).join("")}
      <td><b>${r.union}</b> <span class="note">${pct(r.union, r.n)}%</span></td><td>${r.inter} <span class="note">${pct(r.inter, r.n)}%</span></td>
      <td><a href="#" class="vn-none-lin" data-lin="${escHtml(L)}" style="color:var(--bad);font-weight:600;">${r.none}</a> <span class="note">${pct(r.none, r.n)}%</span></td></tr>`;
  }).join("");
  const tot = [...c.perLin.values()].reduce((a, r) => ({ n: a.n + r.n, union: a.union + r.union, inter: a.inter + r.inter, none: a.none + r.none, each: a.each.map((v, i) => v + r.each[i]) }), { n: 0, union: 0, inter: 0, none: 0, each: c.sel.map(() => 0) });
  const foot = `<tr style="font-weight:700;border-top:2px solid var(--border);"><td>전체</td><td>${tot.n}</td>${tot.each.map((v) => `<td>${v} <span class="note">${pct(v, tot.n)}%</span></td>`).join("")}<td>${tot.union} <span class="note">${pct(tot.union, tot.n)}%</span></td><td>${tot.inter} <span class="note">${pct(tot.inter, tot.n)}%</span></td><td>${tot.none} <span class="note">${pct(tot.none, tot.n)}%</span></td></tr>`;
  return `<div style="overflow:auto;"><table class="mini-reference-table">${head}${body}${foot}</table></div>`;
}

function vnRegionList(c) {
  if (VN.region === null) return "<p class='note'>다이어그램의 숫자나 위 표의 행을 누르면 그 영역의 유전체 목록이 여기에 나옵니다. \"못 잡음\"을 누르면 어떤 세트도 못 잡는 유전체가 나옵니다.</p>";
  const arr = (c.byMask.get(VN.region) || []).slice().sort((a, b) => a.lin.localeCompare(b.lin) || a.id.localeCompare(b.id));
  const byLin = {};
  arr.forEach((r) => { byLin[r.lin] = (byLin[r.lin] || 0) + 1; });
  return `<div class="card" style="margin-top:12px;"><h4 style="margin:0;">${escHtml(vnRegionName(VN.region, c.sel))} <span class="note">· ${arr.length}개${VN.lineage === "ALL" ? "" : ` (${escHtml(VN.lineage)})`}</span></h4>
    <p class="note">${Object.entries(byLin).map(([l, k]) => `${escHtml(l)} ${k}`).join(" · ") || "해당 유전체 없음"}</p>
    <div style="margin:6px 0;"><button type="button" class="btn btn-secondary" id="vnCopyBtn">&#128203; accession 복사</button> <button type="button" class="btn btn-secondary" id="vnCsvBtn">&#128190; CSV 저장</button></div>
    <div style="max-height:320px;overflow:auto;"><table class="mini-reference-table"><tr><th>Accession</th><th>Lineage</th><th>나라</th></tr>
      ${arr.slice(0, 500).map((r) => `<tr><td class="mono">${escHtml(r.id)}</td><td>${escHtml(r.lin)}</td><td>${escHtml(r.country || "—")}</td></tr>`).join("")}</table>
      ${arr.length > 500 ? `<p class="note">처음 500개만 표시합니다. 복사·CSV에는 전부 들어갑니다.</p>` : ""}</div></div>`;
}

function vnDraw(view) {
  const d = VN.data, c = vnCompute();
  view._vn = c;
  const pct = (v) => (c.n ? (100 * v / c.n).toFixed(1) : "0.0");
  const tiles = c.sel.length ? `<div style="display:flex;gap:14px;flex-wrap:wrap;margin:10px 0;">
      <div class="card" style="padding:8px 14px;"><div class="note">대상 유전체</div><b style="font-size:20px;">${c.n}</b></div>
      <div class="card" style="padding:8px 14px;"><div class="note">합집합 (하나라도 잡음)</div><b style="font-size:20px;">${c.union}</b> <span class="note">${pct(c.union)}%</span></div>
      <div class="card" style="padding:8px 14px;"><div class="note">교집합 (모두가 잡음)</div><b style="font-size:20px;">${c.inter}</b> <span class="note">${pct(c.inter)}%</span></div>
      <div class="card" style="padding:8px 14px;cursor:pointer;border-color:var(--bad);" id="vnNoneTile"><div class="note">못 잡음 (여집합)</div><b style="font-size:20px;color:var(--bad);">${c.none}</b> <span class="note">${pct(c.none)}%</span></div></div>` : "";
  const diagram = c.sel.length === 0 ? "<p class='note'>위에서 비교할 프라이머 세트·그룹을 고르세요.</p>"
    : c.sel.length > 4 ? "<p class='note'>5개 이상 고르면 벤 다이어그램 대신 아래 \"겹침 조합 표\"로 보여 줍니다 (4개 이하로 줄이면 그림이 나옵니다).</p>" : vnDiagram(c);
  view.querySelector("#vnBody").innerHTML = `${tiles}
    <div style="display:grid;grid-template-columns:minmax(300px,1.2fr) minmax(280px,1fr);gap:16px;align-items:start;">
      <div class="card vn-diagram-card">${c.sel.length ? vnLegendHtml(c) : ""}${diagram}</div>
      <div class="card"><h4 style="margin:0 0 6px;">겹침 조합 표</h4>${c.sel.length ? vnRegionTable(c) : ""}</div></div>
    <div id="vnList">${vnRegionList(c)}</div>
    <h3 style="margin:18px 0 6px;">리니지별 커버리지 <span class="note">(선택한 세트의 실제 검증 통과 기준 · 리니지 이름을 누르면 위 다이어그램이 그 리니지만 보여 줍니다)</span></h3>
    ${c.sel.length ? vnLineageTable(c) : ""}
    ${vnCurveHtml(c)}`;
}

// ---- compare two Venn combinations ("current selection", a saved one as it was when saved, or a saved one on today's data) ----
// Each side is the UNION of its sets/groups (the panel's coverage); both sides are counted on the genomes they have in common.
// ---- cumulative coverage curves (ROC-style): add the selected sets/groups one by one, best first, and draw how each lineage's coverage rises ----
const VN_LIN_COLORS = ["#4C78A8", "#F58518", "#54A24B", "#E45756", "#B279A2", "#72B7B2", "#9D755D", "#EECA3B", "#FF9DA6", "#79706E"];

function vnCurveData(c) { return vnCurveDataFrom(VN.data.genomes, c.sel); }
function vnCurveDataFrom(genomes, sel) {
  const d = { genomes }, items = sel.map((x) => ({ name: x.name, ids: new Set(x.ids) }));
  const lins = [...new Set(d.genomes.map((g) => g[1]))].sort((a, b) => a.localeCompare(b));
  const total = { ALL: d.genomes.length }, linOf = new Map(d.genomes.map((g) => [g[0], g[1]]));
  lins.forEach((L) => { total[L] = d.genomes.filter((g) => g[1] === L).length; });
  const covered = new Set(), steps = [{ label: "0", gain: 0, cov: Object.fromEntries(["ALL", ...lins].map((k) => [k, 0])) }];
  const left = items.slice();
  while (left.length) {
    let bi = 0, bg = -1;
    left.forEach((it, i) => { let g = 0; it.ids.forEach((id) => { if (!covered.has(id)) g += 1; }); if (g > bg) { bg = g; bi = i; } });
    const [pick] = left.splice(bi, 1);
    pick.ids.forEach((id) => covered.add(id));
    const cov = { ALL: covered.size };
    lins.forEach((L) => { cov[L] = 0; });
    covered.forEach((id) => { const L = linOf.get(id); if (L in cov) cov[L] += 1; });
    steps.push({ label: pick.name, gain: bg, cov });
  }
  return { lins, total, steps };
}

// monotone cubic path through the points (no overshoot between samples), so the curve reads like a ROC curve
function vnSmoothPath(pts) {
  const n = pts.length;
  if (n < 3) return pts.map((p, i) => `${i ? "L" : "M"}${p[0].toFixed(1)},${p[1].toFixed(1)}`).join(" ");
  const dx = [], m = [], t = [];
  for (let i = 0; i < n - 1; i++) { dx[i] = pts[i + 1][0] - pts[i][0]; m[i] = (pts[i + 1][1] - pts[i][1]) / dx[i]; }
  t[0] = m[0]; t[n - 1] = m[n - 2];
  for (let i = 1; i < n - 1; i++) t[i] = m[i - 1] * m[i] <= 0 ? 0 : (m[i - 1] + m[i]) / 2;
  for (let i = 0; i < n - 1; i++) {
    if (m[i] === 0) { t[i] = 0; t[i + 1] = 0; continue; }
    const a = t[i] / m[i], b = t[i + 1] / m[i], s = a * a + b * b;
    if (s > 9) { const k = 3 / Math.sqrt(s); t[i] = k * a * m[i]; t[i + 1] = k * b * m[i]; }
  }
  let path = `M${pts[0][0].toFixed(1)},${pts[0][1].toFixed(1)}`;
  for (let i = 0; i < n - 1; i++) {
    const h = dx[i] / 3;
    path += ` C${(pts[i][0] + h).toFixed(1)},${(pts[i][1] + t[i] * h).toFixed(1)} ${(pts[i + 1][0] - h).toFixed(1)},${(pts[i + 1][1] - t[i + 1] * h).toFixed(1)} ${pts[i + 1][0].toFixed(1)},${pts[i + 1][1].toFixed(1)}`;
  }
  return path;
}

function vnCurveSvg(data) {
  const { lins, total, steps } = data, W = 760, H = 400, L = 54, R = 120, T = 16, B = 96, k = steps.length - 1;
  const X = (i) => L + (k ? (i / k) * (W - L - R) : 0), Y = (p) => T + (1 - p / 100) * (H - T - B);
  const series = [{ key: "ALL", name: "전체", color: "var(--fg)", w: 3.2 }, ...lins.filter((l) => total[l] >= 3).map((l, i) => ({ key: l, name: l, color: VN_LIN_COLORS[i % VN_LIN_COLORS.length], w: 2.2 }))];
  const grid = [0, 25, 50, 75, 100].map((p) => `<line x1="${L}" x2="${W - R}" y1="${Y(p)}" y2="${Y(p)}" stroke="var(--border)" stroke-dasharray="${p % 50 ? "2 4" : "0"}"/><text x="${L - 8}" y="${Y(p) + 4}" text-anchor="end" font-size="11" fill="currentColor" opacity=".75">${p}%</text>`).join("");
  const xl = steps.map((s, i) => `<line x1="${X(i)}" x2="${X(i)}" y1="${H - B}" y2="${H - B + 5}" stroke="currentColor" opacity=".5"/>
    <text transform="translate(${X(i)},${H - B + 14}) rotate(${k > 1 ? 32 : 0})" text-anchor="${k > 1 ? "start" : "middle"}" font-size="11" fill="currentColor">${i ? "+ " : ""}${escHtml(s.label.length > 22 ? s.label.slice(0, 21) + "…" : s.label)}</text>`).join("");
  const diag = `<line x1="${X(0)}" y1="${Y(0)}" x2="${X(k)}" y2="${Y(100)}" stroke="currentColor" opacity=".35" stroke-dasharray="6 5"/>`;
  const lines = series.map((s) => {
    const pts = steps.map((st, i) => [X(i), Y(100 * st.cov[s.key] / (total[s.key] || 1))]);
    const dots = steps.map((st, i) => `<circle cx="${pts[i][0]}" cy="${pts[i][1]}" r="${s.key === "ALL" ? 4.2 : 3.2}" fill="${s.color}" stroke="var(--bg)" stroke-width="1.2"><title>${escHtml(s.name)}: ${i ? `"${escHtml(st.label)}"까지 ${i}개 추가` : "추가 전"} → ${st.cov[s.key]} / ${total[s.key]} (${(100 * st.cov[s.key] / (total[s.key] || 1)).toFixed(1)}%)</title></circle>`).join("");
    return `<path d="${vnSmoothPath(pts)}" fill="none" stroke="${s.color}" stroke-width="${s.w}" stroke-linecap="round"/>${dots}`;
  }).join("");
  const legend = series.map((s) => `<span class="vn-chip"><i style="background:${s.color};${s.key === "ALL" ? "border:1px solid var(--bg);" : ""}"></i>${escHtml(s.name)} <span class="note">n=${total[s.key]}</span></span>`).join("");
  return `<div class="vn-legend">${legend}<span class="vn-chip"><span style="display:inline-block;width:22px;border-top:2px dashed currentColor;opacity:.4;"></span> 같은 속도로 늘 때의 기준선</span></div>
    <svg viewBox="0 0 ${W} ${H}" style="width:100%;max-width:820px;height:auto;display:block;color:var(--fg);">
      <rect x="${L}" y="${T}" width="${W - L - R}" height="${H - T - B}" fill="none" stroke="var(--border)"/>${grid}${diag}${lines}${xl}
      <text x="${(L + W - R) / 2}" y="${H - 6}" text-anchor="middle" font-size="12" fill="currentColor">세트·그룹을 하나씩 더할 때 (새로 가장 많이 잡는 순서)</text>
      <text transform="translate(14,${(T + H - B) / 2}) rotate(-90)" text-anchor="middle" font-size="12" fill="currentColor">리니지 안에서의 커버리지</text></svg>`;
}

function vnCurveHtml(c) {
  if (!c.sel.length) return "";
  const data = vnCurveData(c), steps = data.steps;
  const order = steps.slice(1).map((s, i) => `<li><b>${escHtml(s.label)}</b> <span class="note">+${s.gain}개 새로 잡음 → 전체 ${(100 * s.cov.ALL / (data.total.ALL || 1)).toFixed(1)}%</span></li>`).join("");
  return `<h3 style="margin:18px 0 6px;">리니지별 누적 커버리지 곡선 <span class="note">(ROC 곡선처럼: 위로 빨리 올라갈수록 적은 수의 세트로 많이 잡습니다)</span></h3>
    <div class="card vn-curve-card">${vnCurveSvg(data)}
      <details style="margin-top:8px;"><summary class="note">추가 순서 보기</summary><ol style="margin:6px 0 0 18px;">${order}</ol></details></div>`;
}

function vnCmpSource(val, savedList) {
  const d = VN.data, genomeInfo = new Map(d.genomes.map((g) => [g[0], g]));
  const byKey = new Map(d.items.map((it) => [it.key, it]));
  const live = (keys) => keys.map((k) => byKey.get(k)).filter(Boolean).map((it) => ({ name: it.name, ids: new Set(it.ids) }));
  const everyone = () => new Set(d.genomes.map((g) => g[0]));
  let label, items, pop, lineage = "ALL", note = "";
  if (val === "cur") {
    label = "현재 선택"; items = live([...VN.selected]); pop = everyone(); lineage = VN.lineage;
  } else {
    const [id, mode] = val.split("|"), rec = savedList.find((x) => x.id === id);
    if (!rec) return null;
    lineage = rec.lineage || "ALL";
    if (mode === "snap" && rec.snapshot) {
      label = `${rec.name} (저장 당시 ${rec.created || ""})`;
      items = Object.values(rec.snapshot.items).map((v) => ({ name: v.name, ids: new Set(v.ids) })); pop = new Set(rec.snapshot.pop);
      note = rec.snapshot.scope ? `저장 당시 숙주 기준: ${PD_SCOPE_LABEL[rec.snapshot.scope] || rec.snapshot.scope}` : "";
    } else {
      label = `${rec.name} (현재 데이터)`; items = live(rec.keys); pop = everyone();
    }
  }
  if (lineage !== "ALL") pop = new Set([...pop].filter((id) => (genomeInfo.get(id) || [])[1] === lineage));
  const union = new Set();
  items.forEach((it) => it.ids.forEach((id) => { if (pop.has(id)) union.add(id); }));
  return { label: lineage === "ALL" ? label : `${label} · ${lineage}`, items, pop, union, note };
}

function vnCmpCompute(A, B) {
  const info = new Map(VN.data.genomes.map((g) => [g[0], g]));
  const recs = [...A.pop].filter((id) => B.pop.has(id)).map((id) => ({
    id, lin: (info.get(id) || [])[1] || "NA", country: (info.get(id) || [])[2] || "", m: (A.union.has(id) ? 1 : 0) | (B.union.has(id) ? 2 : 0) }));
  const byMask = new Map();
  recs.forEach((r) => { if (!byMask.has(r.m)) byMask.set(r.m, []); byMask.get(r.m).push(r); });
  const cnt = (m) => (byMask.get(m) || []).length;
  return { sel: [{ name: A.label, ids: [...A.union] }, { name: B.label, ids: [...B.union] }], byMask, n: recs.length, recs, none: cnt(0), full: 3,
    aOnly: cnt(1), bOnly: cnt(2), both: cnt(3), aTotal: cnt(1) + cnt(3), bTotal: cnt(2) + cnt(3), aPop: A.pop.size, bPop: B.pop.size };
}

function vnCmpHtml(c, A, B) {
  const pct = (v, n) => (n ? (100 * v / n).toFixed(1) : "0.0");
  const delta = c.bTotal - c.aTotal, dp = (100 * delta / (c.n || 1));
  const tile = (label, body, extra = "") => `<div class="card" style="padding:8px 14px;${extra}"><div class="note">${label}</div>${body}</div>`;
  const tiles = `<div style="display:flex;gap:14px;flex-wrap:wrap;margin:10px 0;">
    ${tile("공통 대상 유전체", `<b style="font-size:20px;">${c.n}</b> <span class="note">A ${c.aPop} · B ${c.bPop}</span>`)}
    ${tile(`<i class="vn-cmp-sw" style="background:${VN_COLORS[0]}"></i>A 커버리지`, `<b style="font-size:20px;">${c.aTotal}</b> <span class="note">${pct(c.aTotal, c.n)}%</span>`)}
    ${tile(`<i class="vn-cmp-sw" style="background:${VN_COLORS[1]}"></i>B 커버리지`, `<b style="font-size:20px;">${c.bTotal}</b> <span class="note">${pct(c.bTotal, c.n)}%</span>`)}
    ${tile("B − A", `<b style="font-size:20px;color:${delta > 0 ? "var(--ok)" : delta < 0 ? "var(--bad)" : "inherit"};">${delta > 0 ? "+" : ""}${delta}</b> <span class="note">${delta > 0 ? "+" : ""}${dp.toFixed(1)}%p</span>`)}
    ${tile("B만 추가로 잡음", `<b style="font-size:20px;color:var(--ok);">${c.bOnly}</b>`)}
    ${tile("A만 잡음 (B가 놓침)", `<b style="font-size:20px;color:var(--bad);">${c.aOnly}</b>`)}
    ${tile("둘 다 못 잡음 (여집합)", `<b style="font-size:20px;">${c.none}</b> <span class="note">${pct(c.none, c.n)}%</span>`)}</div>`;
  const diagram = vnDiagram(c).replace('id="vnSvg"', 'id="vnCmpSvg"').replace(/vn-region-label/g, "vn-cmp-label");
  const table = vnRegionTable(c).replace(/vn-region-row/g, "vn-cmp-row");
  const lins = [...new Set(c.recs.map((r) => r.lin))].sort((a, b) => a.localeCompare(b));
  const rows = lins.map((L) => {
    const rs = c.recs.filter((r) => r.lin === L), n = rs.length, a = rs.filter((r) => r.m & 1).length, b = rs.filter((r) => r.m & 2).length;
    const ao = rs.filter((r) => r.m === 1).length, bo = rs.filter((r) => r.m === 2).length, no = rs.filter((r) => r.m === 0).length, dd = b - a;
    return `<tr><td>${escHtml(L)}</td><td>${n}</td><td>${a} <span class="note">${pct(a, n)}%</span></td><td>${b} <span class="note">${pct(b, n)}%</span></td>
      <td style="font-weight:700;color:${dd > 0 ? "var(--ok)" : dd < 0 ? "var(--bad)" : "inherit"};">${dd > 0 ? "+" : ""}${dd} <span class="note">(${dd > 0 ? "+" : ""}${(100 * dd / (n || 1)).toFixed(1)}%p)</span></td>
      <td>${ao}</td><td>${bo}</td><td>${no}</td></tr>`;
  }).join("");
  const lineageTable = `<div style="overflow:auto;"><table class="mini-reference-table"><tr><th>${escHtml((PROFILE.lineage && PROFILE.lineage.label) || "Lineage")}</th><th>유전체</th>
    <th style="border-bottom:3px solid ${VN_COLORS[0]};">A 커버리지</th><th style="border-bottom:3px solid ${VN_COLORS[1]};">B 커버리지</th><th>B − A</th><th>A만</th><th>B만</th><th>둘 다 못 잡음</th></tr>${rows}</table></div>`;
  const notes = [A.note && `A: ${A.note}`, B.note && `B: ${B.note}`].filter(Boolean).join(" · ");
  return `<div class="vn-legend"><span class="vn-chip"><i style="background:${VN_COLORS[0]}"></i><b>A</b> ${escHtml(A.label)} <span class="note">(${A.items.length}개 묶음 · 합집합)</span></span>
      <span class="vn-chip"><i style="background:${VN_COLORS[1]}"></i><b>B</b> ${escHtml(B.label)} <span class="note">(${B.items.length}개 묶음 · 합집합)</span></span></div>
    ${notes ? `<p class="note">${escHtml(notes)}</p>` : ""}${c.n < Math.max(c.aPop, c.bPop) ? `<p class="note">두 쪽의 대상 유전체가 달라서(범위·시점 차이) 공통인 ${c.n}개로만 비교했습니다.</p>` : ""}
    ${tiles}
    <div style="display:grid;grid-template-columns:minmax(300px,1.2fr) minmax(280px,1fr);gap:16px;align-items:start;">
      <div class="card">${diagram}</div><div class="card"><h4 style="margin:0 0 6px;">영역별 유전체</h4>${table}</div></div>
    <div id="vnCmpList">${VN.cmpRegion === null || VN.cmpRegion === undefined ? "<p class='note'>다이어그램의 숫자나 표의 행을 누르면 그 영역의 유전체 목록이 나옵니다 (\"B만\" = 새로 잡힌 것, \"A만\" = 놓치게 된 것).</p>" : vnCmpListHtml(c)}</div>
    <h4 style="margin:14px 0 6px;">${escHtml((PROFILE.lineage && PROFILE.lineage.label) || "Lineage")}별 비교</h4>${lineageTable}`;
}

function vnCmpListHtml(c) {
  const arr = (c.byMask.get(VN.cmpRegion) || []).slice().sort((a, b) => a.lin.localeCompare(b.lin) || a.id.localeCompare(b.id));
  return `<div class="card" style="margin:10px 0;"><h4 style="margin:0;">${escHtml(vnRegionName(VN.cmpRegion, c.sel))} <span class="note">· ${arr.length}개</span></h4>
    <div style="margin:6px 0;"><button type="button" class="btn btn-secondary" id="vnCmpCopy">&#128203; accession 복사</button></div>
    <div style="max-height:280px;overflow:auto;"><table class="mini-reference-table"><tr><th>Accession</th><th>Lineage</th><th>나라</th></tr>
    ${arr.slice(0, 400).map((r) => `<tr><td class="mono">${escHtml(r.id)}</td><td>${escHtml(r.lin)}</td><td>${escHtml(r.country || "—")}</td></tr>`).join("")}</table>${arr.length > 400 ? "<p class='note'>처음 400개만 표시합니다. 복사에는 전부 들어갑니다.</p>" : ""}</div></div>`;
}

function vnCmpWire(root, savedList) {
  const selA = root.querySelector("#vnCmpA"), selB = root.querySelector("#vnCmpB"), body = root.querySelector("#vnCmpBody");
  let last = null;
  const draw = () => {
    const A = vnCmpSource(selA.value, savedList), B = vnCmpSource(selB.value, savedList);
    if (!A || !B) { body.innerHTML = "<p class='note'>비교할 쪽을 고르세요.</p>"; return; }
    if (!A.items.length || !B.items.length) { body.innerHTML = "<p class='note'>선택이 비어 있거나 저장된 세트·그룹이 지금은 없어서 비교할 수 없습니다.</p>"; return; }
    last = vnCmpCompute(A, B);
    body.innerHTML = vnCmpHtml(last, A, B);
  };
  selA.addEventListener("change", () => { VN.cmpRegion = null; draw(); });
  selB.addEventListener("change", () => { VN.cmpRegion = null; draw(); });
  body.addEventListener("click", (e) => {
    const reg = e.target.closest(".vn-cmp-label, .vn-cmp-row");
    if (reg && last) { VN.cmpRegion = Number(reg.dataset.mask); draw(); body.querySelector("#vnCmpList").scrollIntoView({ behavior: "smooth", block: "nearest" }); return; }
    if (e.target.closest("#vnCmpCopy") && last) {
      const ids = (last.byMask.get(VN.cmpRegion) || []).map((r) => r.id);
      copyTextToClipboard(ids.join("\n")).then(() => uiAlert(`${ids.length}개 accession을 복사했습니다.`)).catch(() => uiAlert("복사하지 못했습니다."));
    }
  });
  draw();
}

async function showVennPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("vennView");
  const view = $("#vennView");
  const segs = Object.keys((state.hostInfo || {}).segments || {}).map((x) => x.toUpperCase());
  const allSegs = segs.length ? segs : ["S", "L"];
  if (!VN.seg || !allSegs.includes(VN.seg)) VN.seg = allSegs[0];
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try { VN.data = await fetchJSON(`/api/venn/${VN.seg.toLowerCase()}`); } catch (e) { view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`; return; }
  const d = VN.data;
  if (!VN.selected || VN.selSeg !== VN.seg) {                       // first visit / other segment: the groups, else the first sets
    const groups = d.items.filter((it) => it.kind === "group").slice(0, 3);
    VN.selected = new Set((groups.length >= 2 ? groups : d.items.slice(0, 2)).map((it) => it.key));
    VN.selSeg = VN.seg; VN.lineage = "ALL"; VN.region = null;
  }
  const lins = [...new Set(d.genomes.map((g) => g[1]))].sort((a, b) => a.localeCompare(b));
  const byKey = new Map(d.items.map((it) => [it.key, it]));
  const itemChk = (it, level) => {
    const idx = [...VN.selected].indexOf(it.key);
    const col = idx >= 0 ? VN_COLORS[idx % VN_COLORS.length] : "transparent";
    return `<label class="vn-pick${idx >= 0 ? " on" : ""}${it.kind === "group" ? " vn-pick-group" : ""}" style="--c:${col};margin-left:${level > 1 ? (level - 1) * 16 : 0}px;"><input type="checkbox" class="vn-item" value="${escHtml(it.key)}" ${VN.selected.has(it.key) ? "checked" : ""}>
      <span class="vn-dot" style="background:${col};border-color:${idx >= 0 ? col : "var(--muted)"};"></span>
      <span class="vn-name">${it.kind === "group" ? "&#9633; " : ""}${escHtml(it.name)}</span><span class="vn-count">${it.ids.length}</span></label>`;
  };
  // the same group > set > variant tree as the sidebar: pick a whole group, or just its member sets
  const selectorHtml = buildSetTree(VN.seg, "").map((node) => {
    const row = (en, lv) => (byKey.has(en.s.set_id) ? itemChk(byKey.get(en.s.set_id), lv) : "") + en.variants.map((v) => (byKey.has(v.set_id) ? itemChk(byKey.get(v.set_id), lv + 1) : "")).join("");
    if (node.type === "group") return `<div class="vn-tree-group">${byKey.has(node.g.group_id) ? itemChk(byKey.get(node.g.group_id), 0) : ""}<div class="vn-members">${node.kids.map((k) => row(k, 1)).join("")}</div></div>`;
    return `<div class="vn-tree-group vn-tree-loose"><div class="vn-tree-vendor">그룹 없음 · ${escHtml(node.name)}</div><div class="vn-members">${node.sets.map((k) => row(k, 1)).join("")}</div></div>`;
  }).join("") || "<span class='note'>없음</span>";
  let savedList = [];
  try { savedList = (await fetchJSON("/api/venn-saved")).filter((x) => String(x.seg).toUpperCase() === VN.seg.toUpperCase()); } catch (e) { /* none */ }
  const cmpOpts = [`<option value="cur">현재 선택 (${VN.selected.size}개)</option>`].concat(savedList.flatMap((x) => [
    x.snapshot ? `<option value="${x.id}|snap">💾 ${escHtml(x.name)} — 저장 당시 (${escHtml(x.created || "")})</option>` : "",
    `<option value="${x.id}|now">💾 ${escHtml(x.name)} — 현재 데이터</option>`])).join("");
  view.innerHTML = `<h2>🔵 커버리지 벤 다이어그램</h2>
    <p class="note">고른 프라이머 세트·그룹이 <b>실제 검증을 통과시킨 유전체</b>를 겹쳐 봅니다: 합집합(하나라도 잡음), 교집합(모두 잡음), 영역별 단독·겹침, 그리고 어느 것도 못 잡는 유전체(여집합 = 대상 전체 − 합집합). 그룹은 교차 결합으로 통과한 유전체까지 포함합니다. ${pdScopeBadge(d.host_scope && d.host_scope.scope)} 숙주 기준을 따릅니다.</p>
    <div class="seg-tab-row">${allSegs.map((x) => `<button class="seg-tab ${x === VN.seg ? "active" : ""}" data-vnseg="${x}">${x}</button>`).join("")}</div>
    <div class="card">
      <div class="vn-selector-head"><b>비교할 세트·그룹을 고르세요</b> <span class="note">고른 순서대로 색이 정해집니다 · 지금 ${VN.selected.size}개 선택</span></div>
      <div class="vn-selector">${selectorHtml}</div>
      <div class="vn-saved" style="margin-top:8px;"><b class="note">저장한 선택</b> ${savedList.length ? savedList.map((x) => `<span class="vn-saved-chip"><a href="#" class="vn-load" data-id="${x.id}" title="${escHtml(x.keys.length)}개 · ${escHtml(x.created || "")}${x.lineage && x.lineage !== "ALL" ? " · " + escHtml(x.lineage) : ""}">${escHtml(x.name)}</a><button type="button" class="vn-del vn-write" data-id="${x.id}" title="삭제">&times;</button></span>`).join("") : "<span class='note'>아직 없음 — 아래 \"선택 저장\"으로 지금 고른 조합을 이름 붙여 남깁니다</span>"}</div>
      <div class="table-toolbar" style="margin-top:8px;gap:12px;">
        <label class="note">범위 <select id="vnLin"><option value="ALL">전체</option>${lins.map((L) => `<option value="${escHtml(L)}" ${VN.lineage === L ? "selected" : ""}>${escHtml(L)} (${d.genomes.filter((g) => g[1] === L).length})</option>`).join("")}</select></label>
        <button type="button" class="btn btn-secondary" id="vnClear">선택 해제</button>
        <button type="button" class="btn btn-secondary vn-write" id="vnSave" title="지금 고른 세트·그룹과 범위(리니지)를 이름을 붙여 저장합니다 (이 대시보드를 쓰는 모두에게 보입니다)">&#128190; 선택 저장</button>
        <button type="button" class="btn btn-secondary" id="vnPng" title="벤 다이어그램(범례 포함)을 PNG 이미지로 저장합니다">&#128247; 그림 저장</button>
        <span class="note">그림은 4개까지, 더 고르면 표로 보여 줍니다</span></div></div>
    <div id="vnBody"></div>
    <div class="card vn-write" id="vnCompare" style="margin-top:18px;">
      <h3 style="margin:0 0 4px;">&#128256; 벤 다이어그램 비교</h3>
      <p class="note" style="margin:0 0 8px;">두 조합(세트·그룹 묶음)이 잡는 유전체를 비교합니다. 각 조합은 묶음의 <b>합집합</b>(그 패널이 잡는 유전체)입니다. "저장 당시"와 "현재 데이터"를 같은 조합에 골라 두면 <b>재분석 전후의 변화</b>를 볼 수 있습니다.</p>
      <div class="table-toolbar" style="gap:12px;flex-wrap:wrap;">
        <label class="note"><b style="color:${VN_COLORS[0]};">A</b> <select id="vnCmpA">${cmpOpts}</select></label>
        <label class="note"><b style="color:${VN_COLORS[1]};">B</b> <select id="vnCmpB">${cmpOpts}</select></label></div>
      <div id="vnCmpBody"></div></div>`;
  const redraw = () => vnDraw(view);
  view._vnRedraw = redraw;
  view.querySelectorAll("[data-vnseg]").forEach((b) => b.addEventListener("click", () => { VN.seg = b.dataset.vnseg; VN.selected = null; showVennPage(); }));
  view.querySelectorAll(".vn-item").forEach((x) => x.addEventListener("change", () => {
    if (x.checked) VN.selected.add(x.value); else VN.selected.delete(x.value);
    VN.region = null; showVennPage();                                  // re-render so the colour chips follow the selection order
  }));
  view.querySelector("#vnLin").addEventListener("change", (e) => { VN.lineage = e.target.value; VN.region = null; redraw(); });
  view.querySelector("#vnClear").addEventListener("click", () => { VN.selected = new Set(); VN.region = null; showVennPage(); });
  view.querySelector("#vnSave").addEventListener("click", async () => {
    if (!VN.selected.size) { uiAlert("저장할 세트·그룹을 먼저 고르세요."); return; }
    const names = [...VN.selected].map((k) => (byKey.get(k) || {}).name).filter(Boolean);
    const name = await uiPrompt("이 선택의 이름 (같은 이름이 있으면 덮어씁니다)", names.slice(0, 3).join(" + "), { title: "벤 다이어그램 선택 저장" });
    if (!name || !name.trim()) return;
    try {
      await fetchJSON("/api/venn-saved", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ name: name.trim(), seg: VN.seg.toLowerCase(), keys: [...VN.selected], lineage: VN.lineage,
          snapshot: { scope: (d.host_scope || {}).scope || "", pop: d.genomes.map((g) => g[0]),
            items: Object.fromEntries([...VN.selected].map((k) => byKey.get(k)).filter(Boolean).map((it) => [it.key, { name: it.name, ids: it.ids }])) } }) });
      toast("Venn 조합을 저장했습니다.", { type: "ok" });
      showVennPage();
    } catch (e) { uiAlert("저장하지 못했습니다: " + e.message); }
  });
  view.querySelector("#vnPng").addEventListener("click", async () => {
    const card = view.querySelector(".vn-diagram-card");
    if (!card || !VN.selected.size) { uiAlert("그릴 선택이 없습니다."); return; }
    try { downloadBlob(await domToPngBlob(card), `venn_${VN.seg}_${[...VN.selected].join("_")}.png`.slice(0, 120)); } catch (e) { uiAlert(e.message); }
  });
  view.querySelectorAll(".vn-load").forEach((a) => a.addEventListener("click", (e) => {
    e.preventDefault();
    const it = savedList.find((x) => x.id === a.dataset.id);
    if (!it) return;
    VN.selected = new Set(it.keys.filter((k) => byKey.has(k)));
    VN.lineage = it.lineage || "ALL"; VN.region = null;
    showVennPage();
  }));
  view.querySelectorAll(".vn-del").forEach((b) => b.addEventListener("click", async () => {
    if (!(await uiConfirm("이 저장한 선택을 삭제할까요?", { title: "삭제", okLabel: "삭제" }))) return;
    try { await fetchJSON(`/api/venn-saved/${b.dataset.id}/delete`, { method: "POST", headers: { "Content-Type": "application/json" }, body: "{}" }); showVennPage(); } catch (e) { uiAlert(e.message); }
  }));
  {                                                              // comparison panel: A = the newest saved combination as it was (else the current one), B = the current selection
    const root = view.querySelector("#vnCompare"), sA = root.querySelector("#vnCmpA"), sB = root.querySelector("#vnCmpB");
    const withSnap = savedList.filter((x) => x.snapshot);
    sA.value = withSnap.length ? `${withSnap[withSnap.length - 1].id}|snap` : "cur";
    sB.value = withSnap.length ? `${withSnap[withSnap.length - 1].id}|now` : "cur";
    VN.cmpRegion = null;
    vnCmpWire(root, savedList);
  }
  if (!view._vnWired) {                                         // the page element is reused: wire its delegated handlers only once
  view._vnWired = true;
  view.addEventListener("click", (e) => {
    const c = view._vn;
    const reg = e.target.closest(".vn-region-label, .vn-region-row");
    if (reg) { VN.region = Number(reg.dataset.mask); redraw(); view.querySelector("#vnList").scrollIntoView({ behavior: "smooth", block: "nearest" }); return; }
    if (e.target.closest("#vnNoneTile")) { VN.region = 0; redraw(); return; }
    const pl = e.target.closest(".vn-pick-lin");
    if (pl) { e.preventDefault(); VN.lineage = pl.dataset.lin; VN.region = null; view.querySelector("#vnLin").value = VN.lineage; redraw(); return; }
    const nl = e.target.closest(".vn-none-lin");
    if (nl) { e.preventDefault(); VN.lineage = nl.dataset.lin; VN.region = 0; view.querySelector("#vnLin").value = VN.lineage; redraw(); return; }
    const arr = () => (c.byMask.get(VN.region) || []).slice().sort((a, b) => a.lin.localeCompare(b.lin) || a.id.localeCompare(b.id));
    if (e.target.closest("#vnCopyBtn")) { copyTextToClipboard(arr().map((r) => r.id).join("\n")).then(() => uiAlert(`${arr().length}개 accession을 복사했습니다.`)).catch(() => uiAlert("복사하지 못했습니다.")); return; }
    if (e.target.closest("#vnCsvBtn")) {
      const csv = "accession,lineage,country\n" + arr().map((r) => [r.id, r.lin, r.country].map((x) => `"${String(x).replace(/"/g, '""')}"`).join(",")).join("\n") + "\n";
      const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([csv], { type: "text/csv" })); a.download = `venn_${VN.seg}_${vnRegionName(VN.region, c.sel).replace(/[^\w가-힣-]+/g, "_")}.csv`;
      document.body.appendChild(a); a.click(); a.remove(); return;
    }
  });
  // a click on the empty part of a circle (not on a number) also picks the region under the cursor
  view.addEventListener("click", (e) => {
    const svg = e.target.closest("#vnSvg"); if (!svg || e.target.closest(".vn-region-label")) return;
    const c = view._vn; if (!c || !c.sel.length) return;
    const { es, scale } = vnEllipses(c.sel.length), pt = svg.createSVGPoint(); pt.x = e.clientX; pt.y = e.clientY;
    const q = pt.matrixTransform(svg.getScreenCTM().inverse());
    VN.region = vnMask(es, (q.x - 280) / scale, (q.y - 215) / scale); redraw();
  });
  }
  redraw();
}
$("#vennBtn").addEventListener("click", showVennPage);

// ---------------------------------------------------------------------------------------
// 세트 상세 페이지: 프라이머/프로브 결합 위치를 리니지별 pixy π 커브 위에 오버레이 --
// "이 세트의 결합 부위가 population 전체에서 보존된(다양성이 낮은) 구간에 있는가"를
// 한눈에 보여준다. diversityState.data는 리니지별 다양성 페이지와 공유 캐시.
const DIVERSITY_ROLE_META = {
  F: { label: "FORWARD", color: "#4C78A8" },
  P: { label: "PROBE", color: "#F2A65A" },
  R: { label: "REVERSE", color: "#E45756" },
};

async function ensureDiversityData() {
  if (!diversityState.data) {
    diversityState.data = await fetchJSON("/api/diversity");
  }
  return diversityState.data;
}

function bindingOverlaysFromSet(setData) {
  const binding = setData.binding || [];
  return binding
    .filter((b) => b.found && b.start !== null && b.start !== undefined && b.start !== "")
    .map((b) => {
      const meta = DIVERSITY_ROLE_META[b.role] || { label: b.role, color: "#888" };
      // b.set_id is only present on a GROUP's combined binding list (one member set's row
      // tagged with which member it came from) -- a plain set's own binding rows never
      // have it, so the label there stays exactly as before.
      const label = b.set_id ? `${meta.label} (${b.set_name || b.set_id})` : meta.label;
      return { start: Number(b.start), end: Number(b.end), label, color: meta.color, role: b.role, strand: b.strand };
    });
}

function piColor(pi, maxPi) {
  if (pi == null) return "transparent";
  const t = Math.max(0, Math.min(1, pi / maxPi));
  // 3-stop gradient: conserved(low pi, green) -> mid(orange) -> variable(high pi, red) --
  // mirrors the ok/warn/bad palette used everywhere else in this app.
  const stops = [
    [84, 162, 74],   // --ok
    [242, 166, 90],  // --warn/accent2
    [228, 87, 86],   // --bad
  ];
  const seg = t < 0.5 ? 0 : 1;
  const localT = t < 0.5 ? t / 0.5 : (t - 0.5) / 0.5;
  const a = stops[seg], b = stops[seg + 1];
  const r = Math.round(a[0] + (b[0] - a[0]) * localT);
  const g = Math.round(a[1] + (b[1] - a[1]) * localT);
  const bl = Math.round(a[2] + (b[2] - a[2]) * localT);
  return `rgba(${r},${g},${bl},.55)`;
}

function buildAggregatedDiversityWindows(segData) {
  // Averages pi across all real lineages per window -- this baseline row is meant as a quick
  // "is this general area conserved across the population" glance, not a per-lineage
  // breakdown (which is what the dedicated overlay chart above already shows). ALL_LASV is
  // excluded here since it's a pooled-population stat, not another lineage to average in.
  const { real: lineages } = splitLineageKeys(segData);
  if (!lineages.length) return [];
  const base = segData.windows[lineages[0]];
  return base.map((w, i) => {
    let sum = 0, n = 0;
    for (const lin of lineages) {
      const wv = segData.windows[lin][i];
      if (wv && wv.pi != null) { sum += wv.pi; n++; }
    }
    return { start: w.start, end: w.end, pi: n ? sum / n : null };
  });
}

// Column -> reference-genome position for the MSA's pi bars (null for insertion columns that have no reference position).
// The alignment is shown with F on the left, so when F binds the reference's minus strand the columns run DOWN the reference
// coordinates (position = amplicon end - column), not up. Columns where the consensus is a gap are insertions of some genomes
// and do not advance the reference position. The F/Probe/R binding sites (reference coordinates) and their MSA columns anchor an
// exact piecewise-linear map over those ungapped columns. Without usable anchors it falls back to "amplicon start + column".
function msaColumnToRefPos(setData, ampStart) {
  const msa = setData.msa || {};
  const bind = Object.fromEntries((setData.binding || []).filter((b) => b.found !== false).map((b) => [b.role, b]));
  const cols = { F: msa.f_col, P: msa.p_col, R: msa.r_col };
  const fallback = (c) => ampStart + c;
  const strand = bind.F ? bind.F.strand : (bind.R ? (bind.R.strand === "-" ? "+" : "-") : null);
  if (!strand) return fallback;
  const d = strand === "+" ? 1 : -1;
  const cons = msa.consensus || "";
  const u = [];                                           // u[c] = number of non-gap consensus columns before column c
  for (let c = 0, k = 0; c < (msa.columns || 0); c++) { u.push(k); if (cons[c] !== "-") k++; }
  const isGap = (c) => cons[c] === "-";
  const pts = [];
  for (const role of ["F", "P", "R"]) {
    const b = bind[role], cc = cols[role];
    if (!b || !cc) continue;
    const s = Number(b.start), e = Number(b.end), c0 = cc[0], c1 = cc[1] - 1;
    if (!isFinite(s) || !isFinite(e) || !(c1 > c0) || Math.abs((u[c1] - u[c0] + 1) - (e - s + 1)) > 2) continue;   // oligo spans must agree
    pts.push([u[c0], d > 0 ? s : e], [u[c1], d > 0 ? e : s]);
  }
  pts.sort((a, b) => a[0] - b[0]);
  const clean = [];
  for (const p of pts) {                                  // keep only points that run in the expected direction
    const last = clean[clean.length - 1];
    if (!last || (p[0] > last[0] && (p[1] - last[1]) * d > 0)) clean.push(p);
  }
  if (clean.length < 2) return fallback;
  return (c) => {
    if (isGap(c)) return null;
    const x = u[c];
    if (x <= clean[0][0]) return clean[0][1] + d * (x - clean[0][0]);
    const last = clean[clean.length - 1];
    if (x >= last[0]) return last[1] + d * (x - last[0]);
    let i = 1;
    while (clean[i][0] < x) i++;
    const [a0, p0] = clean[i - 1], [a1, p1] = clean[i];
    return p0 + (p1 - p0) * (x - a0) / (a1 - a0);
  };
}

// Ruler text for the MSA header: one character per alignment column. A "|" plus the reference-genome coordinate is written at
// every column whose reference position is a multiple of 10 (insertion columns have none), so the numbers read in the same
// direction as the alignment -- decreasing to the right when F binds the reference's minus strand. Without a position map it is a
// plain "|" every 10 columns.
function buildMsaRuler(columns, posOf) {
  const out = new Array(columns).fill(" ");
  if (!posOf) {
    for (let c = 0; c < columns; c += 10) out[c] = "|";
    return out.join("");
  }
  let free = 0, last = null;
  for (let c = 0; c < columns; c++) {
    const p = posOf(c);
    if (p === null || p === undefined) continue;
    const r = Math.round(p);
    if (r % 10 !== 0 || r === last || c < free) continue;
    const label = "|" + r;
    for (let k = 0; k < label.length && c + k < columns; k++) out[c + k] = label[k];
    free = c + label.length + 1;
    last = r;
  }
  return out.join("");
}

// The pi bars are drawn in an SVG whose width must equal the sequence text above/below it. "N columns x 1ch" is only true when the
// browser's glyph advance equals its `ch` unit; with another font or pixel-rounded glyph advances (e.g. Consolas on Windows) the
// bars came out shorter than the sequences. So the SVG is sized from the real rendered width of a full-length sequence row.
function fitMsaDiversitySvg() {
  const svg = document.querySelector("#msaDiversityRow .msa-diversity-svg");
  const ref = document.querySelector(".msa-consensus .msa-seq") || document.querySelector(".msa-genome-row .msa-seq");
  if (!svg || !ref) return;
  const rng = document.createRange();
  rng.selectNodeContents(ref);
  const w = rng.getBoundingClientRect().width;
  if (w > 0) svg.style.width = `${w}px`;
}
window.addEventListener("resize", fitMsaDiversitySvg);
document.addEventListener("fullscreenchange", () => setTimeout(fitMsaDiversitySvg, 50));

function updateMsaDiversityBaseline(setData, diversityData) {
  const row = $("#msaDiversityRow");
  if (!row) return;
  const seqEl = row.querySelector(".msa-seq");
  const msa = setData.msa;
  const seg = (setData.segment || "").toLowerCase();
  const segData = (diversityData.segments || {})[seg];
  const ampStart = setData.gene && setData.gene.amplicon_start !== undefined && setData.gene.amplicon_start !== null && setData.gene.amplicon_start !== ""
    ? Number(setData.gene.amplicon_start) : null;
  if (!segData || !msa || !msa.columns || ampStart == null) {
    seqEl.innerHTML = "<span class='note'>이 세트에는 다양성 오버레이를 적용할 수 없습니다 (diversity 데이터 또는 결합 위치 정보 없음).</span>";
    return;
  }
  const windowSize = segData.window_size;
  const agg = buildAggregatedDiversityWindows(segData);
  if (!agg.length) {
    seqEl.innerHTML = "<span class='note'>다양성 데이터 없음</span>";
    return;
  }
  const maxPi = Math.max(0.01, ...agg.map((w) => w.pi || 0));

  // Nextstrain의 genome-wide entropy 패널처럼 -- 얇은 색상 띠보다 막대그래프(히스토그램)가
  // 눈에 훨씬 잘 들어온다는 피드백에 따라, window(=column 그룹) 단위로 막대를 그린다.
  // 색 경계가 아니라 window 경계(idx)로 나눠야 인접 window의 값이 우연히 비슷한 색이어도
  // 막대 높이가 뭉개지지 않는다.
  const bars = [];
  const posOf = msaColumnToRefPos(setData, ampStart);
  let curIdx = null, curStart = 0;
  for (let c = 0; c < msa.columns; c++) {
    const pos = posOf(c);
    const idx = pos === null ? "gap" : Math.floor((Math.round(pos) - 1) / windowSize);
    if (idx !== curIdx) {
      if (curIdx !== null) bars.push({ start: curStart, end: c, idx: curIdx });
      curIdx = idx; curStart = c;
    }
  }
  bars.push({ start: curStart, end: msa.columns, idx: curIdx });

  const barH = 100; // viewBox 단위(%) -- preserveAspectRatio="none"으로 실제 픽셀 높이에 맞춰 늘어남
  const barsSvg = bars.map((b) => {
    const w = b.idx === "gap" ? null : (agg[b.idx] || null);
    const pi = w ? w.pi : null;
    const pct = pi != null ? Math.max(3, (pi / maxPi) * 100) : 0;
    const color = piColor(pi, maxPi);
    const piLabel = pi != null ? `π&#8776;${pi.toFixed(3)}` : "N/A";
    const pa = posOf(b.start), pb = posOf(b.end - 1);
    const isGapBar = b.idx === "gap";
    const posLabel = isGapBar ? "" : `~${Math.min(Math.round(pa), Math.round(pb))}-${Math.max(Math.round(pa), Math.round(pb))}`;
    const width = Math.max(0, b.end - b.start);
    if (pi == null) {
      return `<rect x="${b.start}" y="${barH - 1}" width="${width}" height="1" fill="var(--border)"><title>${isGapBar ? "정렬의 삽입 열 (참조 유전체에 없는 위치)" : `위치 ${posLabel} (근사) &#960; 데이터 없음`}</title></rect>`;
    }
    return `<rect x="${b.start}" y="${barH - pct}" width="${width}" height="${pct}" fill="${color}"><title>위치 ${posLabel} (근사) 리니지 평균 ${piLabel}</title></rect>`;
  }).join("");

  seqEl.innerHTML = `<svg viewBox="0 0 ${msa.columns} ${barH}" preserveAspectRatio="none" class="msa-diversity-svg" style="width:${msa.columns}ch;">${barsSvg}</svg>`;
  fitMsaDiversitySvg();
  if (document.fonts && document.fonts.ready) document.fonts.ready.then(fitMsaDiversitySvg);
  const labelEl = row.querySelector(".msa-label");
  if (labelEl) labelEl.innerHTML = `&#960; (다양성)<span class="note" style="font-size:9.5px;display:block;">0&ndash;${maxPi.toFixed(2)}</span>`;
  const legendNote = $("#msaDiversityLegendNote");
  if (legendNote) {
    legendNote.innerHTML = `<span class="dot" style="background:rgba(84,162,74,.55)"></span>&#8594;<span class="dot" style="background:rgba(228,87,86,.55)"></span> &#960;(다양성) 행: 초록=보존(낮은 π) &#8594; 빨강=다양(높은 π), pixy 리니지 평균, window=${windowSize}bp 근사 매핑`;
  }
}

function renderSetDiversitySection(setData) {
  return `
    <div class="card" style="margin-bottom:20px;">
      <h3>리니지별 유전체 다양성(pixy &#960;) 위에 결합 위치 오버레이</h3>
      <div class="note" style="margin-bottom:10px;">
        이 세트의 Forward/Probe/Reverse가 참조유전체(${escHtml(setData.segment)} segment)의 어느 위치에
        결합하는지, 그 위치가 리니지별로 얼마나 보존되어 있는지(&#960;가 낮을수록 보존)를 함께 보여줍니다.
        결합 부위가 다양성이 높은(&#960; 피크) 구간과 겹친다면 일부 리니지/변이체에서 결합이 약해질
        위험이 있다는 뜻입니다.
      </div>
      <div id="setDiversityChartWrap"><p class="note">불러오는 중...</p></div>
    </div>`;
}

function piForRange(segData, start, end) {
  // Average pi per (real) lineage across a reference-coordinate span -- used both for the
  // set-detail overlay's own scaling and to auto-flag a binding site that sits in a locally
  // elevated-diversity patch relative to that lineage's own genome-wide average.
  const { real: lineages } = splitLineageKeys(segData);
  const out = {};
  for (const lin of lineages) {
    const windows = segData.windows[lin] || [];
    let sum = 0, n = 0;
    for (const w of windows) {
      if (w.end < start || w.start > end || w.pi == null) continue;
      sum += w.pi; n++;
    }
    out[lin] = n ? sum / n : null;
  }
  return out;
}

function computePiMismatchLinks(setData, segData, overlays) {
  // Cross-references pixy pi against the genome-by-genome mismatch table this set already
  // computed, so a "this binding site is unusually diverse in lineage X" observation is only
  // surfaced when it's backed by actual observed mismatches in that lineage -- not just a
  // pi bump that mfeprimer's mismatch tolerance happened to absorb without consequence.
  const genomes = setData.genomes || [];
  const genomeWidePi = segData.genome_wide_pi || {};
  const mmField = { F: "f_mismatch", P: "p_mismatch", R: "r_mismatch" };
  const findings = [];
  for (const ov of overlays) {
    const lo = Math.min(ov.start, ov.end), hi = Math.max(ov.start, ov.end);
    const sitePi = piForRange(segData, lo, hi);
    const field = mmField[ov.role];
    for (const [lin, pi] of Object.entries(sitePi)) {
      if (pi == null) continue;
      const baseline = genomeWidePi[lin];
      if (baseline == null || baseline <= 0) continue;
      const ratio = pi / baseline;
      if (pi < 0.03 || ratio < 1.5) continue; // needs to be both a real bump and non-trivial in absolute terms
      const lineageGenomes = genomes.filter((g) => g.lineage === lin);
      // A completely failed match (mismatch field is null -- mfeprimer couldn't align it at
      // all within its own search tolerance) is at least as strong a signal as a small
      // positive mismatch count, so both count as "affected" here.
      const withMismatch = lineageGenomes.filter((g) => g[field] === null || Number(g[field]) > 0);
      if (!lineageGenomes.length || !withMismatch.length) continue;
      findings.push({
        label: ov.label, lineage: lin, pi, baseline, ratio,
        mismatchCount: withMismatch.length, totalCount: lineageGenomes.length,
      });
    }
  }
  return findings;
}

function renderPiMismatchFindings(findings) {
  return findings.map((f) => {
    const pct = Math.round((f.mismatchCount / f.totalCount) * 100);
    return `<li class="warn">&#8505;&#65039; ${escHtml(f.label)} 결합 부위가 <b>${escHtml(f.lineage)}</b> 리니지
      자체 평균보다 ${f.ratio.toFixed(1)}배 높은 다양성(&#960; ${f.pi.toFixed(3)} vs 리니지 평균 ${f.baseline.toFixed(3)})을
      보이는 구간과 겹치며, 실제로 이 리니지 ${f.totalCount}개 게놈 중 ${f.mismatchCount}개(${pct}%)에서 이 부위 mismatch가
      관찰됩니다 (π 오버레이 차트에서 직접 확인 가능).</li>`;
  }).join("");
}

// ---------------------------------------------------------------------------------------
// Fine-tuning suggestion: pixy pi is computed per single base (window_size=1bp), so unlike
// the country/year trend charts, a "shift this oligo a few bp" search actually has real
// resolution to work with. The search itself (+-40bp along the population-pooled ALL_LASV
// pi track) AND the candidate sequence at the new position both run server-side
// (webapp/server.py's compute_shift_suggestions()) rather than here, since building the
// candidate needs to read the population VCF directly -- only the server has cheap
// filesystem access to do that. This just renders whatever setData.shift_suggestions says.
// Deliberately doesn't touch or know about Tm/GC%/secondary structure/distance-to-the-
// other-oligos/amplicon size -- those are exactly the things a human still needs to check
// before actually moving a primer, so this is framed as a starting hint, never an auto-fix.
function renderShiftSuggestions(suggestions) {
  if (!suggestions || !suggestions.length) return "";
  return `<div class="note shift-suggestion-box">
    <b>&#128161; 미세조정(fine-tuning) 제안 (참고용, 자동 적용 아님)</b> — 결합 위치를 그대로 둔 채
    &#960;(다양성)만 기준으로 근처(&plusmn;40bp)를 탐색한 결과입니다. 새 위치의 후보 서열은
    population에서 3% 이상 빈도로 관찰되는 대립형질을 IUPAC 코드로 반영해 구성했습니다. Tm/GC 비율/
    이차구조/다른 올리고와의 거리/실제 합성 가능성(축퇴 코드가 많으면 여러 서열의 혼합물이 되어
    특이도가 떨어질 수 있음)은 전혀 반영되지 않았으니, 실제로 서열을 옮기기 전에 반드시 직접
    재검토하세요.
    <ul style="margin:6px 0 0;padding-left:18px;">
      ${suggestions.map((s) => {
        // Reference-coordinate shift direction only equals the oligo's own 5'->3' reading
        // direction on the "+" strand -- a "-" strand oligo reads toward DECREASING
        // coordinate, so the biological up/downstream label has to flip for it.
        const towardHigherCoord = s.shift > 0;
        const isDownstream = s.strand === "-" ? !towardHigherCoord : towardHigherCoord;
        const label = (DIVERSITY_ROLE_META[s.role] || {}).label || s.role;
        return `<li>${escHtml(label)}를 <b>${Math.abs(s.shift)}bp ${isDownstream ? "하류(3&prime; 방향)" : "상류(5&prime; 방향)"}</b>로 옮기면 이 구간 평균 &#960;가
        ${s.current_pi.toFixed(3)} &rarr; ${s.new_pi.toFixed(3)}로 낮아질 것으로 추정됩니다.
        <div class="mono shift-candidate-seq">${escHtml(s.candidate_seq)}</div></li>`;
      }).join("")}
    </ul>
  </div>`;
}

async function loadSetDiversityOverlay(setData) {
  const wrap = $("#setDiversityChartWrap");
  if (!wrap) return;
  let data;
  try {
    data = await ensureDiversityData();
  } catch (e) {
    wrap.innerHTML = `<p class="note">다양성 데이터를 불러오지 못했습니다: ${e.message}</p>`;
    const row = $("#msaDiversityRow");
    if (row) row.querySelector(".msa-seq").innerHTML = `<span class="note">불러오기 실패: ${escHtml(e.message)}</span>`;
    return;
  }
  updateMsaDiversityBaseline(setData, data);
  const seg = (setData.segment || "").toLowerCase();
  const segData = (data.segments || {})[seg];
  if (!segData) {
    wrap.innerHTML = `<p class="note">아직 ${escHtml(setData.segment)} segment의 diversity 데이터가
      계산되지 않았습니다 (pixy 파이프라인을 실행한 후 다시 확인해주세요).</p>`;
    return;
  }
  const genes = (data.genes || {})[seg] || [];
  const genomeLen = (data.genome_len || {})[seg];
  const overlays = bindingOverlaysFromSet(setData);
  if (!overlays.length) {
    wrap.innerHTML = `<p class="note">참조유전체 결합위치를 찾지 못해 오버레이를 표시할 수 없습니다
      (위 "참조유전체 결합위치" 카드 참고).</p>`;
    return;
  }
  const starts = overlays.map((o) => Math.min(o.start, o.end));
  const ends = overlays.map((o) => Math.max(o.start, o.end));
  const pad = Math.max(40, Math.round((Math.max(...ends) - Math.min(...starts)) * 0.25));
  const xMin = Math.max(1, Math.min(...starts) - pad);
  const xMax = Math.min(genomeLen, Math.max(...ends) + pad);
  const lineages = Object.keys(segData.windows || {}).sort();

  wrap.innerHTML = `
    <div class="diversity-chart-card">
      ${renderDiversityChart(seg, segData, genes, genomeLen, { xMin, xMax, overlays })}
    </div>
    <div class="note" style="margin:6px 0 10px;">결합 부위 주변(&plusmn;${pad}bp)만 확대한 화면입니다.
      전체 genome 구간의 다양성은 <a href="#" id="setDiversityFullLink">리니지별 유전체 다양성</a> 페이지에서
      확인할 수 있습니다.</div>
    ${renderDiversityLegend(lineages, segData.genome_wide_pi || {}, (data.lineage_totals || {})[seg] || {})}
    ${renderShiftSuggestions(setData.shift_suggestions)}`;

  const fullLink = $("#setDiversityFullLink");
  if (fullLink) {
    fullLink.addEventListener("click", (e) => {
      e.preventDefault();
      diversityState.seg = seg;
      showDiversityPage();
    });
  }

  const findings = computePiMismatchLinks(setData, segData, overlays);
  const pointsList = document.querySelector(".assessment-points");
  if (pointsList && findings.length) {
    pointsList.insertAdjacentHTML("beforeend", renderPiMismatchFindings(findings));
  }
}

// ---------------------------------------------------------------------------------------
// 메타데이터 개요 / 지리적 분포 -- population 데이터셋의 출처/수집-필터링 통계와,
// ourworldindata.org/mpox 스타일의 국가별 보고 유전체 수 지도(hover 시 연도별 추이).
const metadataState = { data: null, mapSeg: "total", mapNameToRow: {}, yearChartMode: "line" };

function renderMetadataStatCards(data) {
  const segLabels = { s: "S segment", l: "L segment" };
  const cards = Object.keys(data.segments || {}).map((segKey) => {
    const seg = data.segments[segKey];
    const keptPct = seg.raw_total ? (seg.kept / seg.raw_total) * 100 : 0;
    return `<div class="meta-stat-card">
      <div class="meta-stat-title">${escHtml(segLabels[segKey] || segKey.toUpperCase())}</div>
      <div class="meta-stat-big">${seg.kept.toLocaleString()} <span class="meta-stat-sub">/ ${seg.raw_total.toLocaleString()} (${keptPct.toFixed(0)}%)</span></div>
      <div class="meta-stat-label">QC 통과(필터링 후) / 전체 수집</div>
      <div class="meta-stat-breakdown">
        <div>최소 길이 미달 (&lt;${seg.min_len}bp)<b>${seg.dropped_short.toLocaleString()}</b></div>
        <div>과도한 N/모호염기<b>${seg.dropped_ambiguousN.toLocaleString()}</b></div>
        <div>중복 accession<b>${seg.dropped_dup_accession.toLocaleString()}</b></div>
        <div>메타데이터 없음<b>${seg.dropped_no_metadata.toLocaleString()}</b></div>
      </div>
      <div class="meta-stat-fetched">수집 시점: ${seg.fetched_at ? escHtml(seg.fetched_at.replace("T", " ")) : "알 수 없음"}</div>
    </div>`;
  }).join("");
  return `<div class="meta-stat-row">${cards}</div>`;
}

function renderMetadataCountryTable(data) {
  const rows = data.by_country.slice(0, 20);
  return `<table class="mini-table meta-country-table">
    <thead><tr><th>국가</th><th>S</th><th>L</th><th>합계</th></tr></thead>
    <tbody>${rows.map((r) => `<tr><td>${escHtml(r.country)}</td><td>${r.s}</td><td>${r.l}</td><td><b>${r.total}</b></td></tr>`).join("")}</tbody>
  </table>`;
}

function renderMetadataPageShell(data) {
  const src = data.source || {};
  return `
    <h2>메타데이터 개요 / 지리적 분포</h2>
    <p class="note">population 유전체 데이터의 출처, 수집/필터링 현황과 국가별 보고 현황을 보여줍니다.
      ${src.base_url ? `데이터 출처: <a href="${escHtml(src.base_url)}" target="_blank" rel="noopener">${escHtml(src.base_url)}</a> (Nextstrain)` : ""}</p>
    ${renderMetadataStatCards(data)}
    <h3 style="margin-top:22px;">국가별 보고 유전체 수 (지도)</h3>
    <p class="note" style="margin-top:0;">국가 위에 마우스를 올리면 연도별 보고 추이를 확인할 수 있습니다
      (<a href="https://ourworldindata.org/mpox" target="_blank" rel="noopener">ourworldindata.org/mpox</a> 스타일 참고).</p>
    <div class="seg-tab-row" id="metaMapSegTabs">
      <button class="seg-tab active" data-seg="total">전체</button>
      <button class="seg-tab" data-seg="s">S segment</button>
      <button class="seg-tab" data-seg="l">L segment</button>
    </div>
    <div class="meta-map-card">
      <div id="metaMapSvg" class="meta-map-svg-wrap"></div>
      <div id="metaMapLegend" class="meta-map-legend"></div>
    </div>
    <div id="metaMapTooltip" class="meta-map-tooltip" hidden></div>
    <h3 style="margin-top:22px;">연도별 보고 유전체 수 (추이)</h3>
    <p class="note" style="margin-top:0;">10년 단위로 배경색을 다르게 표시했습니다 (위 지도 탭의 전체/S/L 선택이 함께 적용됩니다).</p>
    <div class="meta-map-card">
      <div id="metaYearChart"></div>
    </div>
    <h3 style="margin-top:22px;">국가별 상세 (상위 20개국)</h3>
    ${renderMetadataCountryTable(data)}
  `;
}

// One distinct background tint per decade (cycled) so 10-year spans are visually separable
// behind the trend lines -- low alpha so lines/points/dark-mode text stay legible on top.
const DECADE_BAND_COLORS = [
  "rgba(66,133,244,.10)", "rgba(52,168,83,.10)", "rgba(251,188,5,.12)",
  "rgba(234,67,53,.10)", "rgba(155,89,182,.10)", "rgba(0,172,193,.10)",
];
const YEAR_CHART_COLOR = d3.scaleOrdinal(d3.schemeTableau10);

// Built with d3 (rather than the hand-rolled SVG the other small trend charts on this page
// use) because it's a multi-series chart with a toggleable legend -- d3's scales/line
// generator/data-join carry that combination better than string-templating would.
// metadataState.yearChartHidden (a Set of series keys) persists which lines are toggled
// off across re-renders (segment-tab switches, legend clicks) for this page visit.
function renderMetaYearChart(data, container) {
  const wrap = document.createElement("div");
  const toggle = document.createElement("div");
  toggle.className = "seg-tab-row";
  toggle.style.marginBottom = "8px";
  toggle.innerHTML = `
    <button class="seg-tab${metadataState.yearChartMode !== "percent" ? " active" : ""}" data-mode="line">절대 건수</button>
    <button class="seg-tab${metadataState.yearChartMode === "percent" ? " active" : ""}" data-mode="percent">리니지 비율 (100% 누적)</button>
  `;
  const chartDiv = document.createElement("div");
  wrap.appendChild(toggle);
  wrap.appendChild(chartDiv);
  container.innerHTML = "";
  container.appendChild(wrap);
  toggle.querySelectorAll(".seg-tab").forEach((btn) => {
    btn.addEventListener("click", () => {
      metadataState.yearChartMode = btn.dataset.mode;
      renderMetaYearChart(data, container);
    });
  });

  if (metadataState.yearChartMode === "percent") {
    renderMetaYearPercentAreaChart(data, chartDiv);
  } else {
    renderMetaYearLineChart(data, chartDiv);
  }
}

function renderMetaYearLineChart(data, container) {
  const seg = metadataState.mapSeg || "total";
  const byYear = data.by_year || {};
  const byLineageYear = data.by_lineage_year || {};
  metadataState.yearChartHidden = metadataState.yearChartHidden || new Set();
  const valueOf = (counts) => (seg === "total" ? (counts.s || 0) + (counts.l || 0) : (counts[seg] || 0));

  const yearSet = new Set(Object.keys(byYear));
  Object.values(byLineageYear).forEach((years) => Object.keys(years).forEach((y) => yearSet.add(y)));
  const years = Array.from(yearSet).map(Number).sort((a, b) => a - b);
  container.innerHTML = "";
  if (!years.length) { container.innerHTML = "<p class='note'>연도 정보 없음</p>"; return; }

  const lineageTotal = (lin) => Object.values(byLineageYear[lin] || {}).reduce((sum, c) => sum + valueOf(c), 0);
  const lineages = Object.keys(byLineageYear).sort((a, b) => lineageTotal(b) - lineageTotal(a));

  const series = [
    { key: "__total__", label: "전체", color: "var(--accent)", bold: true,
      values: years.map((y) => valueOf(byYear[y] || {})) },
    ...lineages.map((lin) => ({
      key: lin, label: lin, color: YEAR_CHART_COLOR(lin), bold: false,
      values: years.map((y) => valueOf((byLineageYear[lin] || {})[y] || {})),
    })),
  ];
  const visible = series.filter((s) => !metadataState.yearChartHidden.has(s.key));
  const drawOrder = [...visible].sort((a, b) => (a.bold ? 1 : 0) - (b.bold ? 1 : 0)); // bold("전체") drawn last, on top

  const w = 900, h = 300, mL = 46, mR = 16, mT = 14, mB = 34;
  const minYear = years[0], maxYear = years[years.length - 1];
  const x = d3.scaleLinear().domain([minYear, maxYear]).range([mL, w - mR]);
  const maxV = Math.max(1, ...(visible.length ? visible.flatMap((s) => s.values) : [0]));
  const y = d3.scaleLinear().domain([0, maxV]).nice().range([h - mB, mT]);

  const svg = d3.select(container).append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`).attr("class", "meta-year-trend-svg");

  const decadeStart = Math.floor(minYear / 10) * 10;
  let colorIdx = 0;
  for (let d = decadeStart; d <= maxYear; d += 10) {
    const bandStart = Math.max(d, minYear), bandEnd = Math.min(d + 10, maxYear + 1);
    if (bandEnd > bandStart) {
      svg.append("rect")
        .attr("x", x(bandStart)).attr("y", mT)
        .attr("width", Math.max(0, x(bandEnd) - x(bandStart))).attr("height", h - mT - mB)
        .attr("fill", DECADE_BAND_COLORS[colorIdx % DECADE_BAND_COLORS.length]);
      svg.append("text")
        .attr("x", (x(bandStart) + x(bandEnd)) / 2).attr("y", h - mB + 14)
        .attr("text-anchor", "middle").attr("font-size", 9).attr("fill", "var(--muted)")
        .text(`${d}s`);
    }
    colorIdx++;
  }

  svg.append("line").attr("x1", mL).attr("y1", h - mB).attr("x2", w - mR).attr("y2", h - mB).attr("stroke", "var(--border)");
  svg.append("line").attr("x1", mL).attr("y1", mT).attr("x2", mL).attr("y2", h - mB).attr("stroke", "var(--border)");
  svg.append("text").attr("x", mL - 6).attr("y", mT + 8).attr("text-anchor", "end").attr("font-size", 9).attr("fill", "var(--muted)").text(maxV);
  svg.append("text").attr("x", mL - 6).attr("y", h - mB).attr("text-anchor", "end").attr("font-size", 9).attr("fill", "var(--muted)").text("0");

  const lineGen = d3.line().x((v, i) => x(years[i])).y((v) => y(v));
  const seriesG = svg.selectAll("g.series").data(drawOrder, (d) => d.key).join("g").attr("class", "series");
  seriesG.append("path")
    .attr("fill", "none").attr("stroke", (d) => d.color)
    .attr("stroke-width", (d) => (d.bold ? 2.6 : 1.6))
    .attr("d", (d) => lineGen(d.values));
  seriesG.each(function pointsForSeries(d) {
    d3.select(this).selectAll("circle").data(d.values).join("circle")
      .attr("cx", (v, i) => x(years[i])).attr("cy", (v) => y(v)).attr("r", d.bold ? 2.8 : 2)
      .attr("fill", d.color)
      .append("title").text((v, i) => `${d.label} ${years[i]}: ${v}개`);
  });

  const legend = document.createElement("div");
  legend.className = "meta-year-legend";
  legend.innerHTML = series.map((s) => {
    const off = metadataState.yearChartHidden.has(s.key);
    return `<span class="meta-year-legend-chip${off ? " off" : ""}" data-key="${escHtml(s.key)}" style="--chip-color:${s.color}">
      <span class="dot"></span>${escHtml(s.label)}
    </span>`;
  }).join("");
  container.appendChild(legend);
  legend.querySelectorAll(".meta-year-legend-chip").forEach((chip) => {
    chip.addEventListener("click", () => {
      const key = chip.dataset.key;
      if (metadataState.yearChartHidden.has(key)) metadataState.yearChartHidden.delete(key);
      else metadataState.yearChartHidden.add(key);
      renderMetaYearLineChart(data, container);
    });
  });
}

// 100%-stacked area view (Highcharts "area-stacked-percent" style) of the same per-year,
// per-lineage data the line chart above uses -- shows how the COMPOSITION (relative share)
// of reported genomes across lineages shifted over time, rather than absolute counts.
// d3.stack().offset(d3.stackOffsetExpand) does the normalize-to-100%-per-year math; a
// lineage toggled off in the shared legend is simply left out of the stack, so the
// remaining ones renormalize to fill 100% again automatically.
function renderMetaYearPercentAreaChart(data, container) {
  const seg = metadataState.mapSeg || "total";
  const byLineageYear = data.by_lineage_year || {};
  metadataState.yearChartHidden = metadataState.yearChartHidden || new Set();
  const valueOf = (counts) => (seg === "total" ? (counts.s || 0) + (counts.l || 0) : (counts[seg] || 0));

  const yearSet = new Set();
  Object.values(byLineageYear).forEach((years) => Object.keys(years).forEach((y) => yearSet.add(y)));
  const years = Array.from(yearSet).map(Number).sort((a, b) => a - b);
  container.innerHTML = "";
  if (!years.length) { container.innerHTML = "<p class='note'>연도 정보 없음</p>"; return; }

  const lineageTotal = (lin) => Object.values(byLineageYear[lin] || {}).reduce((sum, c) => sum + valueOf(c), 0);
  const allLineages = Object.keys(byLineageYear).sort((a, b) => lineageTotal(b) - lineageTotal(a));
  const lineages = allLineages.filter((lin) => !metadataState.yearChartHidden.has(lin));

  // one row per year, one column per visible lineage -- what d3.stack() expects
  const rows = years.map((y) => {
    const row = { year: y };
    lineages.forEach((lin) => { row[lin] = valueOf((byLineageYear[lin] || {})[y] || {}); });
    return row;
  });
  const stack = d3.stack().keys(lineages).order(d3.stackOrderNone).offset(d3.stackOffsetExpand);
  const series = lineages.length ? stack(rows) : [];

  const w = 900, h = 300, mL = 46, mR = 16, mT = 14, mB = 34;
  const minYear = years[0], maxYear = years[years.length - 1];
  const x = d3.scaleLinear().domain([minYear, maxYear]).range([mL, w - mR]);
  const y = d3.scaleLinear().domain([0, 1]).range([h - mB, mT]);

  const svg = d3.select(container).append("svg")
    .attr("viewBox", `0 0 ${w} ${h}`).attr("class", "meta-year-trend-svg");

  const decadeStart = Math.floor(minYear / 10) * 10;
  let colorIdx = 0;
  for (let d = decadeStart; d <= maxYear; d += 10) {
    const bandStart = Math.max(d, minYear), bandEnd = Math.min(d + 10, maxYear + 1);
    if (bandEnd > bandStart) {
      svg.append("rect")
        .attr("x", x(bandStart)).attr("y", mT)
        .attr("width", Math.max(0, x(bandEnd) - x(bandStart))).attr("height", h - mT - mB)
        .attr("fill", DECADE_BAND_COLORS[colorIdx % DECADE_BAND_COLORS.length]);
      svg.append("text")
        .attr("x", (x(bandStart) + x(bandEnd)) / 2).attr("y", h - mB + 14)
        .attr("text-anchor", "middle").attr("font-size", 9).attr("fill", "var(--muted)")
        .text(`${d}s`);
    }
    colorIdx++;
  }

  svg.append("line").attr("x1", mL).attr("y1", h - mB).attr("x2", w - mR).attr("y2", h - mB).attr("stroke", "var(--border)");
  svg.append("line").attr("x1", mL).attr("y1", mT).attr("x2", mL).attr("y2", h - mB).attr("stroke", "var(--border)");
  [0, 0.25, 0.5, 0.75, 1].forEach((t) => {
    svg.append("text").attr("x", mL - 6).attr("y", y(t) + 3).attr("text-anchor", "end").attr("font-size", 9).attr("fill", "var(--muted)").text(`${Math.round(t * 100)}%`);
  });

  const areaGen = d3.area()
    .x((d, i) => x(years[i]))
    .y0((d) => y(d[0])).y1((d) => y(d[1]));
  svg.selectAll("path.area").data(series, (d) => d.key).join("path")
    .attr("class", "area")
    .attr("fill", (d) => YEAR_CHART_COLOR(d.key)).attr("opacity", 0.85)
    .attr("d", areaGen)
    .append("title").text((d) => d.key);

  const legend = document.createElement("div");
  legend.className = "meta-year-legend";
  legend.innerHTML = allLineages.map((lin) => {
    const off = metadataState.yearChartHidden.has(lin);
    return `<span class="meta-year-legend-chip${off ? " off" : ""}" data-key="${escHtml(lin)}" style="--chip-color:${YEAR_CHART_COLOR(lin)}">
      <span class="dot"></span>${escHtml(lin)}
    </span>`;
  }).join("");
  container.appendChild(legend);
  legend.querySelectorAll(".meta-year-legend-chip").forEach((chip) => {
    chip.addEventListener("click", () => {
      const key = chip.dataset.key;
      if (metadataState.yearChartHidden.has(key)) metadataState.yearChartHidden.delete(key);
      else metadataState.yearChartHidden.add(key);
      renderMetaYearPercentAreaChart(data, container);
    });
  });
}

function renderCountryTrendSVG(yearCounts, seg) {
  const years = Object.keys(yearCounts).sort();
  if (!years.length) return "<div class='note' style='margin:4px 0 0;'>연도 정보 없음</div>";
  const w = 220, h = 90, mL = 30, mR = 8, mT = 8, mB = 18;
  const values = years.map((y) => (seg === "total" ? (yearCounts[y].s || 0) + (yearCounts[y].l || 0) : (yearCounts[y][seg] || 0)));
  const maxV = Math.max(...values, 1);
  const xScale = (i) => mL + (years.length > 1 ? (i / (years.length - 1)) * (w - mL - mR) : (w - mL - mR) / 2);
  const yScale = (v) => (h - mB) - (v / maxV) * (h - mT - mB);
  let d = "";
  values.forEach((v, i) => { d += (i === 0 ? "M " : "L ") + xScale(i).toFixed(1) + "," + yScale(v).toFixed(1) + " "; });
  return `<svg viewBox="0 0 ${w} ${h}" width="${w}" height="${h}" class="meta-trend-svg">
    <line x1="${mL}" y1="${h - mB}" x2="${w - mR}" y2="${h - mB}" stroke="var(--border)"></line>
    <path d="${d}" fill="none" stroke="var(--accent)" stroke-width="1.6"></path>
    ${values.map((v, i) => `<circle cx="${xScale(i).toFixed(1)}" cy="${yScale(v).toFixed(1)}" r="2" fill="var(--accent)"></circle>`).join("")}
    <text x="${mL}" y="${h - 4}" font-size="9" fill="var(--muted)">${escHtml(years[0])}</text>
    <text x="${w - mR}" y="${h - 4}" text-anchor="end" font-size="9" fill="var(--muted)">${escHtml(years[years.length - 1])}</text>
    <text x="${mL - 4}" y="${mT + 6}" text-anchor="end" font-size="9" fill="var(--muted)">${maxV}</text>
  </svg>`;
}

function showMetaMapTooltip(event, mapName, data, seg) {
  const row = metadataState.mapNameToRow[mapName];
  if (!row) return;
  const tooltip = $("#metaMapTooltip");
  const yearCounts = (data.by_country_year || {})[row.country] || {};
  const segLabel = seg === "total" ? "전체" : seg.toUpperCase();
  const count = seg === "total" ? row.total : (row[seg] || 0);
  tooltip.innerHTML = `<div class="meta-tooltip-head">${escHtml(mapName)}</div>
    <div class="meta-tooltip-count">${escHtml(segLabel)}: ${count.toLocaleString()}개 게놈</div>
    ${renderCountryTrendSVG(yearCounts, seg)}`;
  tooltip.hidden = false;
  const pad = 14;
  let left = event.clientX + pad, top = event.clientY + pad;
  if (left + 240 > window.innerWidth) left = event.clientX - 240 - pad;
  if (top + 140 > window.innerHeight) top = event.clientY - 140 - pad;
  tooltip.style.left = `${left}px`;
  tooltip.style.top = `${top}px`;
}

function renderMetaMapLegend(maxCount) {
  const legend = $("#metaMapLegend");
  if (!legend || typeof d3 === "undefined") return;
  const stops = d3.range(0, 1.0001, 0.1).map((t) => d3.interpolateYlOrRd(t));
  legend.innerHTML = `<div class="meta-legend-bar" style="background:linear-gradient(to right, ${stops.join(",")})"></div>
    <div class="meta-legend-labels"><span>0</span><span>보고된 유전체 수</span><span>${maxCount.toLocaleString()}</span></div>`;
}

async function renderMetadataMap(data, container) {
  if (typeof d3 === "undefined" || typeof topojson === "undefined") {
    container.innerHTML = "<p class='note'>지도 라이브러리(d3/topojson)를 불러오지 못했습니다.</p>";
    return;
  }
  container.innerHTML = "";
  const width = 900, height = 440;
  const seg = metadataState.mapSeg || "total";

  const mapNameToRow = {};
  const countByMapName = {};
  for (const row of data.by_country) {
    if (!row.map_name) continue;
    mapNameToRow[row.map_name] = row;
    const v = seg === "total" ? row.total : (row[seg] || 0);
    if (v > 0) countByMapName[row.map_name] = v;
  }
  metadataState.mapNameToRow = mapNameToRow;

  let topoData;
  try {
    topoData = await fetchJSON("/vendor/countries-110m.json");
  } catch (e) {
    container.innerHTML = `<p class="note">지도 데이터를 불러오지 못했습니다: ${e.message}</p>`;
    return;
  }

  const countries = topojson.feature(topoData, topoData.objects.countries).features;
  const projection = d3.geoNaturalEarth1().fitSize([width, height], { type: "Sphere" });
  const path = d3.geoPath(projection);
  const maxCount = Math.max(1, ...Object.values(countByMapName));
  const colorScale = d3.scaleSequential(d3.interpolateYlOrRd).domain([0, Math.log1p(maxCount)]);

  const svg = d3.select(container).append("svg")
    .attr("viewBox", `0 0 ${width} ${height}`)
    .attr("class", "meta-map-svg");

  svg.append("path").datum({ type: "Sphere" }).attr("d", path)
    .attr("fill", "var(--card)").attr("stroke", "var(--border)");

  svg.selectAll("path.country")
    .data(countries)
    .join("path")
    .attr("class", "country")
    .attr("d", path)
    .attr("fill", (d) => {
      const c = countByMapName[d.properties.name];
      return c ? colorScale(Math.log1p(c)) : "var(--card)";
    })
    .attr("stroke", "var(--border)")
    .attr("stroke-width", 0.5)
    .style("cursor", (d) => (countByMapName[d.properties.name] ? "pointer" : "default"))
    .on("mousemove", (event, d) => {
      if (!countByMapName[d.properties.name]) { $("#metaMapTooltip").hidden = true; return; }
      showMetaMapTooltip(event, d.properties.name, data, seg);
    })
    .on("mouseleave", () => { $("#metaMapTooltip").hidden = true; });

  renderMetaMapLegend(maxCount);
}

function wireMetadataMapTabs(view) {
  view.querySelectorAll("#metaMapSegTabs .seg-tab").forEach((btn) => {
    btn.addEventListener("click", async () => {
      view.querySelectorAll("#metaMapSegTabs .seg-tab").forEach((b) => b.classList.toggle("active", b === btn));
      metadataState.mapSeg = btn.dataset.seg;
      await renderMetadataMap(metadataState.data, view.querySelector("#metaMapSvg"));
      renderMetaYearChart(metadataState.data, view.querySelector("#metaYearChart"));
    });
  });
}

async function showMetadataPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("metadataView");
  const view = $("#metadataView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    metadataState.data = await fetchJSON("/api/metadata");
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
    return;
  }
  view.innerHTML = renderMetadataPageShell(metadataState.data);
  wireMetadataMapTabs(view);
  await renderMetadataMap(metadataState.data, view.querySelector("#metaMapSvg"));
  renderMetaYearChart(metadataState.data, view.querySelector("#metaYearChart"));
}

$("#metadataBtn").addEventListener("click", showMetadataPage);

async function pollJob(jobId) {
  const jt = jobToast("세트 분석");
  while (true) {
    const job = await fetchJSON(`/api/jobs/${jobId}`);
    if (job.status === "running") {
      jt.progress(job);
      $("#fStatus").textContent = job.progress || "분석 중...";
      await new Promise((r) => setTimeout(r, 1200));
      continue;
    }
    if (job.status === "error") {
      $("#fStatus").textContent = `오류: ${job.error}`;
      jt.done(false, "세트 분석 실패: " + (job.error || ""));
      return null;
    }
    jt.done(true);
    return job.set_id;
  }
}

async function submitNewSet() {
  const payload = {
    name: $("#fName").value.trim(),
    segment: $("#fSegment").value,
    f_seq: $("#fF").value.trim(),
    r_seq: $("#fR").value.trim(),
    p_seq: $("#fP").value.trim(),
    save_to_manifest: $("#fSave").checked,
  };
  if (!payload.name || !payload.f_seq || !payload.r_seq || !payload.p_seq) {
    $("#fStatus").textContent = "이름, F/R/Probe 서열을 모두 입력하세요.";
    return;
  }
  $("#fSubmit").disabled = true;
  $("#fStatus").textContent = "제출 중...";
  try {
    const { job_id, error } = await fetchJSON("/api/analyze", {
      method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload),
    });
    if (error) throw new Error(error);
    const setId = await pollJob(job_id);
    if (setId) {
      await loadSets();
      selectSet(setId);
    }
  } catch (e) {
    $("#fStatus").textContent = `오류: ${e.message}`;
  } finally {
    $("#fSubmit").disabled = false;
  }
}

// ---------------- overview (전체 대시보드) ----------------

function heatmapCellStyle(v, vmax) {
  if (!v) return "";
  const alpha = 0.08 + 0.82 * (v / vmax);
  return `background: rgba(76,120,168,${alpha.toFixed(2)})`;
}

// The overview lists groups and sets in the same tree order as the sidebar (group, then its sets, variants right under their parent,
// then the sets without a group by vendor) instead of by validated %, so a set is always found in the same place.
function overviewOrder(data) {
  const setById = new Map((data.sets || []).map((x) => [x.set_id, x]));
  const groupById = new Map((data.groups || []).map((x) => [x.group_id, x]));
  const rows = [];
  const segs = [...new Set([...(data.sets || []).map((x) => x.segment), ...(data.groups || []).map((x) => x.segment)])].sort((x, y) => (x === "S" ? 0 : x === "L" ? 1 : 2) - (y === "S" ? 0 : y === "L" ? 1 : 2) || x.localeCompare(y));
  const addSet = (en, level) => {
    const s = setById.get(en.s.set_id);
    if (!s) return;
    rows.push({ kind: "set", s, level });
    for (const v of en.variants) { const vs = setById.get(v.set_id); if (vs) rows.push({ kind: "set", s: vs, level: level + 1 }); }
  };
  for (const seg of segs) {
    for (const node of buildSetTree(seg, "")) {
      if (node.type === "group") {
        const g = groupById.get(node.g.group_id);
        if (!g) continue;
        rows.push({ kind: "group", g });
        node.kids.forEach((k) => addSet(k, 1));
      } else {
        rows.push({ kind: "vendor", name: node.name, seg });
        node.sets.forEach((k) => addSet(k, 0));
      }
    }
  }
  return rows;
}

function renderOverviewHtml(data) {
  const order = overviewOrder(data);
  const sets = order.filter((r) => r.kind === "set").map((r) => r.s);
  const levelOf = new Map(order.filter((r) => r.kind === "set").map((r) => [r.s.set_id, r.level]));
  const groups = data.groups || [];
  const lineageRows = data.lineage || [];
  const ds = data.dataset_sizes || {};

  const setRowHtml = (s) => {
    const sm = s.summary || {};
    const lv = levelOf.get(s.set_id) || 0;
    return `<tr class="genome-row" data-set-id="${s.set_id}">
      <td class="mono" style="padding-left:${10 + lv * 16}px;">${lv ? "↳ " : ""}${escHtml(s.display_name || s.set_id)}${s.duplicate_of ? ` <span class="dup-badge" title="F/R/Probe가 ${escHtml(s.duplicate_of)}와 같은 서열">= ${escHtml(s.duplicate_of)}</span>` : ""}</td>
      <td><span class="pill pill-${s.segment}">${s.segment}</span></td>
      <td>${s.source}</td>
      <td>${fmtPct(sm.primer_coverage_mm0_pct)}</td>
      <td>${fmtPct(sm.primer_coverage_mm1_pct)}</td>
      <td>${fmtPct(sm.primer_coverage_mm2_pct)}</td>
      <td>${fmtPct(sm.probe_coverage_mm0_pct)}</td>
      <td>${sm.probe_dropout ?? "—"}</td>
      <td style="font-weight:700;color:var(--accent);">${fmtPct(sm.validated_pct)}</td>
      <td>${evaluationPillHtml(s.evaluation)}</td>
      <td>${finalChipHtml(s.set_id)}</td>
    </tr>`;
  };
  const tableRows = order.map((r) => {
    if (r.kind === "set") return setRowHtml(r.s);
    if (r.kind === "group") return `<tr class="genome-row ov-group-head" data-group-id="${r.g.group_id}"><td colspan="11">&#9633; <b>${escHtml(r.g.label)}</b> <span class="note">(${(r.g.members || "").split(";").filter(Boolean).length}개 세트 ALL · 합집합 ${fmtPct(r.g.validated_pct)} · 교차 포함 ${fmtPct(r.g.validated_with_cross_pct ?? r.g.validated_pct)}) ${finalChipHtml(r.g.group_id)}</span></td></tr>`;
    return `<tr class="ov-vendor-head"><td colspan="11">그룹 없음 · ${escHtml(r.name)}</td></tr>`;
  }).join("");

  const groupRows = groups.map((g) => `
    <tr class="genome-row" data-group-id="${g.group_id}">
      <td class="mono">${g.label}</td>
      <td><span class="pill pill-${g.segment}">${g.segment}</span></td>
      <td style="font-size:11px;">${escHtml((g.member_names || []).join(", "))}</td>
      <td>${g.total_genomes}</td>
      <td style="font-weight:700;color:var(--accent);">${g.validated} (${fmtPct(g.validated_pct)})</td>
      <td style="font-weight:700;color:var(--accent);">${g.validated_with_cross ?? g.validated} (${fmtPct(g.validated_with_cross_pct ?? g.validated_pct)})</td>
      <td>${finalChipHtml(g.group_id)}</td>
    </tr>`).join("");

  // lineage heatmap (count and % of that lineage's total population captured), merged into
  // the SAME row/table as each set's coverage bars -- one glance covers both instead of
  // relating two separately-titled blocks by eye.
  const lineageTotals = data.lineage_totals || {};
  const totals = {};
  const perSet = {};
  for (const r of lineageRows) {
    perSet[r.set_id] = perSet[r.set_id] || {};
    perSet[r.set_id][r.lineage] = Number(r.validated_count);
    totals[r.lineage] = (totals[r.lineage] || 0) + Number(r.validated_count);
  }
  const lineages = Object.keys(totals).sort((a, b) => totals[b] - totals[a]);
  const vmax = Math.max(1, ...lineageRows.map((r) => Number(r.validated_count)));

  // Groups ("ALL" set-combination union coverage) get their own rows in the same merged
  // table -- they only have ONE coverage metric (the union's validated_pct, no separate
  // primer/probe split the way an individual set does), so their coverage cell shows a
  // single bar instead of the usual three.
  const groupHeatRow = (g) => {
    const segTotals = lineageTotals[g.segment] || {};
    const cells = lineages.map((l) => {
      const v = (perSet[g.group_id] || {})[l] || 0;
      const denom = segTotals[l] || 0;
      const pct = denom ? Math.round((v / denom) * 100) : null;
      const label = v ? `${v}${pct !== null ? `<br><span class="hm-pct">${pct}%</span>` : ""}` : "";
      return `<td style="${heatmapCellStyle(v, vmax)}">${label}</td>`;
    }).join("");
    const memberCount = (g.members || "").split(";").filter(Boolean).length;
    return `
      <tr class="cov-heat-row cov-heat-row-group" data-group-id="${g.group_id}">
        <th><span class="pill pill-${g.segment}">${g.segment}</span> &#9633; ${escHtml(g.label)} <span class="hm-total">(${memberCount}개 세트 ALL)</span></th>
        <td class="cov-cell">
          <div class="cov-bars">
            <div class="cov-bar-track" title="세트별 합집합"><div class="cov-bar-fill validated" style="width:${Number(g.validated_pct) || 0}%"></div></div>
            <div class="cov-bar-track" title="교차 결합 포함"><div class="cov-bar-fill probe" style="width:${Number(g.validated_with_cross_pct ?? g.validated_pct) || 0}%"></div></div>
          </div>
        </td>
        ${cells}
      </tr>`;
  };

  const mergedRow = (s) => {
    const sm = s.summary || {};
    const std = s.standard || {};
    const oxHtml = std.ok ? `<span class="ox ox-o">O</span>` : `<span class="ox ox-x">X</span>`;
    const segTotals = lineageTotals[s.segment] || {};
    const cells = lineages.map((l) => {
      const v = (perSet[s.set_id] || {})[l] || 0;
      const denom = segTotals[l] || 0;
      const pct = denom ? Math.round((v / denom) * 100) : null;
      const label = v ? `${v}${pct !== null ? `<br><span class="hm-pct">${pct}%</span>` : ""}` : "";
      return `<td style="${heatmapCellStyle(v, vmax)}">${label}</td>`;
    }).join("");
    return `
      <tr class="cov-heat-row" data-set-id="${s.set_id}">
        <th style="padding-left:${10 + (levelOf.get(s.set_id) || 0) * 16}px;"><span class="pill pill-${s.segment}">${s.segment}</span> ${levelOf.get(s.set_id) ? "↳ " : ""}${escHtml(s.display_name || s.set_id)} ${oxHtml}</th>
        <td class="cov-cell">
          <div class="cov-bars">
            <div class="cov-bar-track"><div class="cov-bar-fill primer" style="width:${sm.primer_coverage_mm1_pct || 0}%"></div></div>
            <div class="cov-bar-track"><div class="cov-bar-fill probe" style="width:${sm.probe_coverage_mm0_pct || 0}%"></div></div>
            <div class="cov-bar-track"><div class="cov-bar-fill validated" style="width:${sm.validated_pct || 0}%"></div></div>
          </div>
        </td>
        ${cells}
      </tr>`;
  };
  const heatBody = order.map((r) => {
    if (r.kind === "group") return groupHeatRow(r.g);
    if (r.kind === "set") return mergedRow(r.s);
    return `<tr class="ov-vendor-head"><td colspan="${2 + lineages.length}">그룹 없음 · ${escHtml(r.name)}</td></tr>`;
  }).join("");
  const lineageHeaderRow = lineages.map((l) => {
    // find denominator from whichever segment this lineage belongs to (S/L totals rarely overlap by name)
    const denom = (lineageTotals.S && lineageTotals.S[l]) ?? (lineageTotals.L && lineageTotals.L[l]) ?? null;
    return `<th>${escHtml(l)}${denom ? `<div class="hm-total">n=${denom}</div>` : ""}</th>`;
  }).join("");

  return `
    <h2>${escHtml(PROFILE.short_name || "LASV")} Primer/Probe 전체 대시보드</h2>
    ${ictvCrumbHtml()}
    <div class="subtitle">Nextstrain(lassa/s, /l) 최신 데이터셋 &middot; 판정 기준: primer mismatch&le;1, probe mismatch=0(완전일치)</div>

    <div class="stat-tiles">
      <div class="stat-tile"><div class="label">S segment 게놈</div><div class="value">${ds.S ?? "—"}</div></div>
      <div class="stat-tile"><div class="label">L segment 게놈</div><div class="value">${ds.L ?? "—"}</div></div>
      <div class="stat-tile"><div class="label">등록된 세트</div><div class="value">${sets.length}</div></div>
      <div class="stat-tile clickable-card" id="stdTile" role="button" tabindex="0" title="클릭하면 세트별 표준물질 매치를 창으로 보여줍니다"><div class="label">표준물질 O</div><div class="value">${sets.some((s) => s.standard?.available !== false) ? `${sets.filter((s) => s.standard?.ok).length} / ${sets.filter((s) => s.standard?.available !== false).length}` : "&mdash;"}</div></div>
    </div>

    <div class="section-title">세트별 커버리지 &amp; 리니지 분포 (primer vs probe vs 최종 validated, n=해당 리니지 전체 게놈 수)</div>
    <div class="cov-legend">
      <span><span class="dot" style="background:#B0B7C1"></span>Primer match (mm&le;1)</span>
      <span><span class="dot" style="background:var(--accent2)"></span>Probe 완전일치</span>
      <span><span class="dot" style="background:var(--accent)"></span>Validated (최종)</span>
    </div>
    <div class="heatmap-wrap">
      <table class="heatmap cov-heatmap"><thead><tr><th></th><th>커버리지</th>${lineageHeaderRow}</tr></thead>
      <tbody>${heatBody}</tbody></table>
    </div>

    <div class="section-title">세트별 상세 요약 (그룹 → 세트 순서, 사이드바와 같은 순서 · 클릭하면 상세보기)</div>
    <div class="note" style="margin:2px 0 8px;">
      <b>평가</b> 기준(각 세트 상세 페이지의 "종합 평가" 박스와 동일한 자동 채점): Population
      커버리지(&ge;95% +2점 / &ge;80% +1점 / &ge;50% -1점 / 미만 -2점), Probe dropout 비율(30%
      초과 -1점, 매우 낮으면 +1점), 표준물질(KDCA) 매치(O +1점 / X -2점), 참조유전체 결합 실패 여부(실패
      시 -2점, 정상이고 amplicon 크기가 적절하면 +1점), 이차구조(dimer/hairpin) 위험(있으면 -1점, 없으면
      +1점)을 합산해 총점 3점 이상이면 <span class="pill pill-ok">통과</span>, 0~2점이면
      <span class="pill pill-warn">보류</span>, 0점 미만이면 <span class="pill pill-bad">실패</span>로
      표시합니다. 참고용 자동 요약이며 최종 판단은 실험 검증을 거쳐야 합니다.
    </div>
    <div class="genome-table-wrap" style="max-height:none;">
      <table class="genome-table">
        <thead><tr>
          <th>프라이머 세트</th><th>segment</th><th>source</th>
          <th>primer mm=0</th><th>primer mm=1</th><th>primer mm=2</th>
          <th>probe mm=0</th><th>probe dropout</th><th>validated</th>
          <th title="세트 상세 페이지의 종합 평가와 동일한 기준 — 위 설명 참고">평가</th>
          <th title="사용자가 세트 상세 화면에서 직접 표시한 최종 결과">최종 결과</th>
        </tr></thead>
        <tbody>${tableRows}</tbody>
      </table>
    </div>

    ${groups.length ? `
    <div class="section-title">세트 조합("ALL") 합집합 커버리지</div>
    <div class="genome-table-wrap" style="max-height:none;">
      <table class="genome-table">
        <thead><tr><th>조합</th><th>segment</th><th>구성 세트</th><th>total</th><th title="세트별 단독 통과의 합집합">합집합 validated</th><th title="합집합 + 세트 간 교차 조합으로 통과하는 게놈">교차 포함 validated</th><th title="사용자가 그룹 상세 화면에서 직접 표시한 최종 결과">최종 결과</th></tr></thead>
        <tbody>${groupRows}</tbody>
      </table>
    </div>` : ""}
  `;
}

async function showOverview() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("overview");
  const overview = $("#overview");
  overview.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    await loadFinalResults();
    const data = await fetchJSON("/api/overview");
    overview.innerHTML = hostScopeBadgeHtml(data.host_scope) + renderOverviewHtml(data);
    overview.querySelectorAll("[data-set-id]").forEach((el2) => {
      el2.addEventListener("click", () => selectSet(el2.dataset.setId));
    });
    overview.querySelectorAll("[data-group-id]").forEach((el2) => {
      el2.addEventListener("click", () => selectGroup(el2.dataset.groupId));
    });
  } catch (e) {
    overview.innerHTML = `<p class="note">불러오기 실패: ${e.message}</p>`;
  }
}

// ---------------- pipeline control (데이터 업데이트 / 분석 다시하기) ----------------

function setPipelineButtonsDisabled(disabled) {
  $("#fetchDataBtn").disabled = disabled;
  $("#rerunBtn").disabled = disabled;
}

function fmtElapsed(sec) {
  if (sec == null) return "";
  const m = Math.floor(sec / 60);
  const s = sec % 60;
  return m > 0 ? `${m}분 ${s}초` : `${s}초`;
}

function showPipelineModal(label) {
  $("#pipelineModalTitle").textContent = label;
  $("#pipelineModalSub").textContent = "시작 중...";
  $("#pipelineModalLog").textContent = "";
  $("#pipelineModal").hidden = false;
}

function updatePipelineModal(label, sub, logText) {
  if ($("#pipelineModal").hidden) return; // user closed it -- job keeps running server-side, just stop touching the DOM
  $("#pipelineModalTitle").textContent = label;
  $("#pipelineModalSub").textContent = sub;
  const logEl = $("#pipelineModalLog");
  const wasAtBottom = logEl.scrollHeight - logEl.scrollTop - logEl.clientHeight < 40;
  logEl.textContent = logText || "(아직 로그 없음)";
  if (wasAtBottom) logEl.scrollTop = logEl.scrollHeight;
}

$("#pipelineModalClose").addEventListener("click", () => { $("#pipelineModal").hidden = true; });

async function pollPipelineJob(jobId, label) {
  // Resolves to the job's final status ("done" / "error") so callers can chain a follow-up step.
  const statusEl = $("#pipelineStatus");
  statusEl.hidden = false;
  let finalStatus = "error";
  const jt = jobToast(label);
  while (true) {
    let job;
    try {
      job = await fetchJSON(`/api/jobs/${jobId}`);
    } catch (e) {
      statusEl.className = "pipeline-status error";
      statusEl.textContent = `${label}: 상태 확인 실패 (${e.message})`;
      jt.done(false, `${label}: 상태 확인 실패`);
      updatePipelineModal(label, `상태 확인 실패: ${e.message}`, $("#pipelineModalLog").textContent);
      break;
    }
    if (job.status === "running") {
      statusEl.className = "pipeline-status running";
      statusEl.textContent = `${label} 실행 중... ${job.progress || ""} (${fmtElapsed(job.elapsed_sec)} 경과)`;
      jt.progress(job);
      updatePipelineModal(label, `실행 중... ${job.progress || ""} (${fmtElapsed(job.elapsed_sec)} 경과)`, job.log_tail);
      await new Promise((r) => setTimeout(r, 1500));
      continue;
    }
    finalStatus = job.status;
    if (job.status === "error") {
      statusEl.className = "pipeline-status error";
      statusEl.textContent = `${label} 실패 (returncode=${job.returncode}). 로그: ${job.log_path || ""}`;
      jt.done(false);
      updatePipelineModal(label, `실패 (returncode=${job.returncode})`, job.log_tail);
    } else {
      statusEl.className = "pipeline-status done";
      statusEl.textContent = `${label} 완료.`;
      jt.done(true);
      updatePipelineModal(label, "완료", job.log_tail);
      await loadSets();
      if (state.activeId) {
        selectSet(state.activeId);
      } else {
        showOverview();
      }
    }
    break;
  }
  setPipelineButtonsDisabled(false);
  return finalStatus;
}

async function startPipelineJob(endpoint, label) {
  setPipelineButtonsDisabled(true);
  const statusEl = $("#pipelineStatus");
  statusEl.hidden = false;
  statusEl.className = "pipeline-status running";
  statusEl.textContent = `${label} 시작 중...`;
  showPipelineModal(label);
  try {
    const res = await fetch(endpoint, { method: "POST" });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    return await pollPipelineJob(data.job_id, label);
  } catch (e) {
    statusEl.className = "pipeline-status error";
    statusEl.textContent = `${label} 시작 실패: ${e.message}`;
    updatePipelineModal(label, `시작 실패: ${e.message}`, "");
    setPipelineButtonsDisabled(false);
    return "error";
  }
}

// ---------------- startup check: is there newer Nextstrain data? ----------------
// Runs lazily after the first render (never blocks loading). If the remote files changed since
// the last 데이터 업데이트, offer to update and then re-run the analysis. "이 버전 건너뛰기"
// remembers the remote fingerprint so the same update isn't asked about on every reload.
async function checkFreshnessOnStartup() {
  let fresh;
  try { fresh = await fetchJSON("/api/pipeline/freshness"); } catch { return; }
  const changed = (fresh.changed_files || []).filter((f) => f.reason === "content_changed");
  if (!changed.length) return;   // "no_baseline" alone = nothing to compare with, don't nag
  const fingerprint = changed.map((f) => `${f.label}:${f.record_count}`).join("|");
  try { if (localStorage.getItem("lasv_fresh_skip") === fingerprint) return; } catch { /* private mode */ }
  try {
    const active = await fetchJSON("/api/jobs/active");
    if (Array.isArray(active) && active.length) return;   // a run is in progress elsewhere -- don't pile on
  } catch { /* fine -- just ask */ }
  const rows = changed.map((f) =>
    `<li>${escHtml(f.label)}: ${f.prev_record_count ?? "?"}개 → <b>${f.record_count ?? "?"}개</b></li>`).join("");
  const choice = await uiDialog({
    title: "새 Nextstrain 데이터가 있습니다",
    html: `<div>마지막 데이터 업데이트 이후 Nextstrain의 원본 파일이 바뀌었습니다.</div><ul class="ui-dialog-list">${rows}</ul>
      <div>데이터를 업데이트한 뒤 새 데이터로 전체 분석을 다시 돌릴까요? (업데이트는 몇 분, 재분석은 더 오래 걸릴 수 있습니다)</div>`,
    buttons: [
      { label: "나중에", cancel: true, value: () => "later" },
      { label: "이 버전 건너뛰기", value: () => "skip" },
      { label: "업데이트 후 재분석", primary: true, value: () => "go" },
    ],
  });
  if (choice === "skip") {
    try { localStorage.setItem("lasv_fresh_skip", fingerprint); } catch { /* private mode */ }
    return;
  }
  if (choice !== "go") return;
  const status = await startPipelineJob("/api/pipeline/fetch-data", "데이터 업데이트");
  if (status !== "done") return;
  if (await uiConfirm("데이터 업데이트가 끝났습니다. 새 데이터로 전체 분석을 다시 실행할까요?", { okLabel: "재분석 시작", title: "재분석" })) {
    await startPipelineJob("/api/pipeline/rerun-analysis", "분석 다시하기");
  }
}

// ---------------- report PDF generation ----------------

async function pollReportJob(jobId) {
  const jt = jobToast("종합 리포트 PDF");
  while (true) {
    let job;
    try {
      job = await fetchJSON(`/api/jobs/${jobId}`);
    } catch (e) {
      updatePipelineModal("종합 리포트 PDF 생성", `상태 확인 실패: ${e.message}`, "");
      break;
    }
    if (job.status === "running") {
      updatePipelineModal("종합 리포트 PDF 생성", `생성 중... ${job.progress || ""} (${fmtElapsed(job.elapsed_sec)} 경과)`, job.log_tail);
      jt.progress(job);
      await new Promise((r) => setTimeout(r, 1200));
      continue;
    }
    if (job.status === "error") {
      updatePipelineModal("종합 리포트 PDF 생성", `실패 (returncode=${job.returncode})`, job.log_tail);
      jt.done(false);
    } else {
      updatePipelineModal("종합 리포트 PDF 생성", "완료!", job.log_tail);
      jt.done(true);
      const sub = $("#pipelineModalSub");
      if (sub && !sub.querySelector("a")) {
        sub.appendChild(document.createTextNode(" "));
        const link = document.createElement("a");
        link.href = "/api/report/download";
        link.textContent = "📄 PDF 다운로드";
        link.className = "btn btn-primary";
        link.style.marginLeft = "8px";
        sub.appendChild(link);
      }
    }
    break;
  }
  $("#reportBtn").disabled = false;
}

async function startReportJob() {
  $("#reportBtn").disabled = true;
  showPipelineModal("종합 리포트 PDF 생성");
  try {
    const res = await fetch("/api/report/generate", { method: "POST" });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    pollReportJob(data.job_id);
  } catch (e) {
    updatePipelineModal("종합 리포트 PDF 생성", `시작 실패: ${e.message}`, "");
    $("#reportBtn").disabled = false;
  }
}

$("#reportBtn").addEventListener("click", startReportJob);

// ---------------- init ----------------

$("#addSetBtn").addEventListener("click", showAddForm);
$("#overviewBtn").addEventListener("click", showOverview);
$("#searchBox").addEventListener("input", renderSidebar);
$("#fSubmit").addEventListener("click", submitNewSet);
$("#fetchDataBtn").addEventListener("click", async () => {
  // Nextstrain의 원본 파일이 바뀌었는지 먼저 HEAD로 가볍게 확인 -- 안 바뀌었으면 굳이 몇 분씩
  // 걸리는 재다운로드+재인덱싱을 할 이유가 없으므로, 그 경우엔 "정말 진행할지" 한 번 더 확인.
  // 확인 자체가 실패하면(네트워크 문제 등) 안전하게 기존 동작(그냥 확인 후 진행)으로 넘어감.
  let confirmMsg = "Nextstrain에서 최신 데이터를 새로 받고 mfeprimer 인덱스를 다시 만듭니다. 계속할까요?";
  try {
    const fresh = await fetchJSON("/api/pipeline/freshness");
    if (!fresh.changed && (fresh.unchanged_files || []).length) {
      const detail = fresh.unchanged_files.map((f) => `${f.label} ${f.record_count || "?"}개`).join(", ");
      confirmMsg = `현재 로컬 데이터가 Nextstrain 최신 버전과 동일한 것으로 보입니다 (변경 없음: ${detail}). `
        + "그래도 다시 다운로드하시겠습니까? (보통은 필요하지 않습니다)";
    }
  } catch (e) {
    // freshness 확인 자체가 실패해도 업데이트 버튼은 계속 동작해야 하므로 조용히 무시.
  }
  if (await uiConfirm(confirmMsg)) {
    startPipelineJob("/api/pipeline/fetch-data", "데이터 업데이트");
  }
});
$("#rerunBtn").addEventListener("click", async () => {
  if (await uiConfirm("현재 데이터로 전체 세트(기존+사용자 추가분 포함) 분석을 다시 계산합니다. 몇 분 이상 걸릴 수 있습니다. 계속할까요?")) {
    startPipelineJob("/api/pipeline/rerun-analysis", "분석 다시하기");
  }
});

function setGroupMode(on) {
  state.groupMode = on;
  if (!on) state.selectedForGroup.clear();
  $("#groupBar").hidden = !on;
  $("#groupModeBtn").textContent = on ? "그룹 만들기 취소" : "+ 그룹 만들기 (여러 세트 묶기)";
  if (on) $("#groupLabelInput").value = "";
  updateGroupBar();
  renderSidebar();
}

$("#groupModeBtn").addEventListener("click", () => setGroupMode(!state.groupMode));
$("#groupCancelBtn").addEventListener("click", () => setGroupMode(false));
$("#groupLabelInput").addEventListener("input", updateGroupBar);
$("#groupCreateBtn").addEventListener("click", async () => {
  const label = $("#groupLabelInput").value.trim();
  const memberIds = [...state.selectedForGroup];
  if (!label || memberIds.length < 2) return;
  $("#groupCreateBtn").disabled = true;
  try {
    const res = await fetch("/api/groups", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ label, member_set_ids: memberIds }),
    });
    const data = await res.json();
    if (!res.ok) throw new Error(data.error || `HTTP ${res.status}`);
    setGroupMode(false);
    await loadGroups();
    if (data.job_id) {
      pollPipelineJob(data.job_id, `그룹 '${label}' 계산`).then(() => loadGroups());
    }
    if (data.warning) uiAlert(data.warning);
  } catch (e) {
    uiAlert(`그룹 생성 실패: ${e.message}`);
    $("#groupCreateBtn").disabled = false;
  }
});

// ---------------- active-jobs banner (visible from any connected computer, not just the
// one that started the job -- the server has no concept of "users", so a second person
// looking at the dashboard from a different machine has no other way to know a
// 데이터 업데이트/분석 다시하기/새 세트 분석/그룹 계산 is already running elsewhere) ----------------

function fmtElapsedShort(sec) {
  if (sec == null) return "";
  const m = Math.floor(sec / 60), s = sec % 60;
  return m > 0 ? `${m}분 ${s}초` : `${s}초`;
}

async function pollActiveJobsBanner() {
  const banner = $("#activeJobsBanner");
  try {
    const jobs = await fetchJSON("/api/jobs/active");
    if (jobs.length) {
      banner.hidden = false;
      banner.innerHTML = jobs.map((j) =>
        `<div class="active-job-row">⏳ ${j.label || (j.kind === "analyze" ? "새 세트 분석" : "파이프라인 작업")} 실행 중… (${fmtElapsedShort(j.elapsed_sec)} 경과)${j.progress ? `<div class="active-job-progress">${j.progress}</div>` : ""}</div>`
      ).join("");
    } else {
      banner.hidden = true;
    }
  } catch {
    // transient network hiccup -- leave whatever was last shown, try again next tick
  }
}
setInterval(pollActiveJobsBanner, 4000);
pollActiveJobsBanner();

// ---------------------------------------------------------------------------------------
// Manual light/dark toggle -- the app already followed the OS's prefers-color-scheme; this
// adds an explicit override (persisted in localStorage, applied pre-paint by a small inline
// script in index.html's <head> so reloading never flashes the wrong theme for a frame).
const THEME_KEY = "lasv_theme";
function currentTheme() {
  return document.documentElement.dataset.theme
    || (matchMedia("(prefers-color-scheme: dark)").matches ? "dark" : "light");
}
function applyTheme(theme) {
  document.documentElement.dataset.theme = theme;
  try { localStorage.setItem(THEME_KEY, theme); } catch { /* private mode etc. */ }
  $("#themeToggleBtn").innerHTML = theme === "dark" ? "&#9728;" : "&#9789;";
}
$("#themeToggleBtn").addEventListener("click", () => applyTheme(currentTheme() === "dark" ? "light" : "dark"));
applyTheme(currentTheme());

// Vaadin-Lumo skin (lumo.css, scoped to html[data-skin="lumo"]) vs the original look.
function applySkin(skin) {
  document.documentElement.dataset.skin = skin;
  try { localStorage.setItem("lasv_skin", skin); } catch { /* private mode etc. */ }
  $("#skinToggleBtn").textContent = skin === "lumo" ? "스타일: Lumo" : "스타일: 기본";
}
$("#skinToggleBtn").addEventListener("click", () => applySkin(document.documentElement.dataset.skin === "lumo" ? "classic" : "lumo"));
applySkin(document.documentElement.dataset.skin === "classic" ? "classic" : "lumo");

$("#langToggleBtn").addEventListener("click", () => setLang(currentLang() === "en" ? "ko" : "en"));
setLang(currentLang());

(async () => {
  await loadProfile();
  await loadUiState();
  await loadHostScope();
  await Promise.all([loadSets(), loadGroups()]);
  showOverview();
  // first visit: offer the UI tour; the data-freshness question waits until that is answered
  setTimeout(async () => { await offerTourOnFirstVisit(); checkFreshnessOnStartup(); }, 1500);   // lazy: after the first paint, never blocks it
})();

// ---------------------------------------------------------------------------------------------------------------------------
// 📚 문헌 프라이머 (sub-menus LASV / Nipah): primer/probe sets collected from papers. The data is config/literature/<virus>.tsv on THIS
// computer only (git-ignored); an entry counts as "원문 대조 완료" only when someone wrote down where it was checked (webapp/literature.py).
const LIT_FIELDS = [
  ["name", "이름"], ["segment", "분절"], ["forward_primer", "Forward (5'→3')"], ["reverse_primer", "Reverse (5'→3')"], ["probe", "Probe (5'→3')"],
  ["probe_dyes", "형광/소광체"], ["target", "표적 유전자·영역"], ["citation", "출처 (저자, 저널, 연도)"], ["doi", "DOI"], ["pmid", "PMID"], ["license", "이용 조건 (예: CC BY 4.0)"],
  ["verified_note", "원문 확인 위치 (예: Table 2)"], ["remark", "메모"],
];
let litState = { payload: null, virus: "lasv", editing: null, errors: [] };

async function showLiteraturePage() {
  litState.editing = null; litState.errors = [];
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("litView");
  const view = $("#litView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  try {
    litState.payload = await fetchJSON("/api/literature");
    litState.virus = litState.payload.edition;
    renderLiteraturePage();
  } catch (e) {
    view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`;
  }
}

function litLinks(r) {
  const out = [];
  if (r.doi) out.push(`<a href="https://doi.org/${escHtml(r.doi)}" target="_blank" rel="noopener noreferrer">DOI</a>`);
  if (r.pmid) out.push(`<a href="https://pubmed.ncbi.nlm.nih.gov/${escHtml(r.pmid)}/" target="_blank" rel="noopener noreferrer">PubMed</a>`);
  return out.join(" · ");
}

function renderLiteraturePage() {
  const view = $("#litView"), p = litState.payload, v = p.viruses[litState.virus];
  const rows = v.rows, mine = true;
  const badge = (r) => r.verified === "yes"
    ? `<span class="lit-badge lit-ok" title="${escHtml(r.verified_note)}">원문 대조 완료</span>`
    : `<span class="lit-badge lit-no">미검증</span>`;
  const seq = (s) => (s ? `<span class="mono lit-seq">${escHtml(s)}</span>` : `<span class="muted">없음</span>`);
  const body = rows.map((r, i) => `<tr>
      <td>${badge(r)}</td>
      <td><b>${escHtml(r.name)}</b><div class="note">${escHtml(r.lit_id)}</div></td>
      <td>${escHtml(r.segment)}</td>
      <td>${seq(r.forward_primer)}</td><td>${seq(r.reverse_primer)}</td><td>${seq(r.probe)}</td>
      <td>${escHtml(r.probe_dyes)}</td><td>${escHtml(r.target)}</td>
      <td>${escHtml(r.citation)}<div class="note">${litLinks(r)}${r.license ? ` · 이용 조건: ${escHtml(r.license)}` : ""}</div>${r.verified_note ? `<div class="note">확인: ${escHtml(r.verified_note)}</div>` : ""}</td>
      <td class="lit-actions">
        ${r.in_analysis ? `<span class="note">분석 목록에 있음</span>`
          : r.can_add ? `<button class="btn btn-primary" data-lit-add="${i}">분석에 추가</button>`
          : `<span class="note" title="${mine ? "Forward/Reverse/Probe 서열이 모두 있어야 합니다" : "다른 바이러스의 세트입니다"}">${mine ? "서열 불완전" : "이 에디션에서 분석 불가"}</span>`}
        <button class="btn btn-secondary" data-lit-edit="${i}">수정</button>
        <button class="btn btn-secondary" data-lit-del="${i}">삭제</button>
      </td></tr>`).join("");
  view.innerHTML = `
    <h2>📚 문헌 프라이머</h2>
    <p class="note">논문에서 찾은 TaqMan 프라이머·프로브 목록입니다. 서열은 <code>${escHtml(v.file)}</code>(이 컴퓨터, git 제외)에만 저장되며 GitHub에 올라가지 않습니다.
    <b>‘원문 대조 완료’는 논문 원문과 서열을 직접 비교하고 확인한 위치를 적은 항목에만 붙습니다.</b> 그 외 서열은 틀렸을 수 있으니 분석에 쓰기 전에 원문을 확인하세요.</p>
    <p><button class="btn btn-primary" id="litNewBtn">+ 항목 추가</button>
    <span class="note" id="litStatus"></span></p>
    <div id="litForm"></div>
    ${rows.length ? `<div class="table-wrap"><table class="primers-table lit-table"><thead><tr>
      <th>상태</th><th>이름</th><th>분절</th><th>Forward</th><th>Reverse</th><th>Probe</th><th>염료</th><th>표적</th><th>출처</th><th></th></tr></thead><tbody>${body}</tbody></table></div>`
      : `<p class="note">아직 항목이 없습니다. 논문에서 확인한 세트를 ‘+ 항목 추가’로 입력하세요.</p>`}`;
  $("#litNewBtn").onclick = () => litEdit(null);
  view.querySelectorAll("[data-lit-edit]").forEach((b) => (b.onclick = () => litEdit(Number(b.dataset.litEdit))));
  view.querySelectorAll("[data-lit-del]").forEach((b) => (b.onclick = () => litDelete(Number(b.dataset.litDel))));
  view.querySelectorAll("[data-lit-add]").forEach((b) => (b.onclick = () => litAddToAnalysis(Number(b.dataset.litAdd))));
  if (litState.editing !== undefined && litState.editing !== false) litRenderForm();
}

function litEdit(i) { litState.editing = i === null ? { index: null, row: {} } : { index: i, row: { ...litState.payload.viruses[litState.virus].rows[i] } }; litState.errors = []; litRenderForm(); }

function litRenderForm() {
  const box = $("#litForm"), e = litState.editing;
  if (!e) { box.innerHTML = ""; return; }
  const err = (f) => litState.errors.filter((x) => x.field === f).map((x) => `<div class="cell-error-text">${escHtml(x.message)}</div>`).join("");
  const field = ([f, label]) => {
    const seqCls = ["forward_primer", "reverse_primer", "probe"].includes(f) ? " mono" : "";
    const val = escHtml(e.row[f] || "");
    const input = f === "segment"
      ? `<select data-lit="segment">${litState.payload.viruses[litState.virus].segments.map((s) => `<option ${e.row.segment === s ? "selected" : ""}>${s}</option>`).join("")}</select>`
      : `<input type="text" class="cell-input${seqCls}" data-lit="${f}" value="${val}" spellcheck="false" autocomplete="off">`;
    return `<label class="lit-field"><span>${label}</span>${input}${err(f)}</label>`;
  };
  box.innerHTML = `<div class="lit-form">
    <h3>${e.index === null ? "새 항목" : "항목 수정"}</h3>
    <div class="lit-grid">${LIT_FIELDS.map(field).join("")}</div>
    <label class="lit-check"><input type="checkbox" data-lit="verified" ${e.row.verified === "yes" ? "checked" : ""}>
      원문과 서열을 직접 비교했습니다 (체크하면 ‘원문 확인 위치’가 필요합니다)</label>
    ${err("verified_note")}
    <p><button class="btn btn-primary" id="litSaveBtn">저장</button> <button class="btn btn-secondary" id="litCancelBtn">취소</button></p></div>`;
  $("#litCancelBtn").onclick = () => { litState.editing = null; litState.errors = []; litRenderForm(); };
  $("#litSaveBtn").onclick = litSaveForm;
  box.scrollIntoView({ block: "nearest" });
}

async function litPostRows(rows) {
  const res = await fetch("/api/literature/save", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ virus: litState.virus, rows }) });
  const data = await res.json();
  if (!res.ok) { litState.errors = data.errors || [{ field: "name", message: data.error || "저장 실패" }]; return false; }
  litState.payload = data.payload; return true;
}

async function litSaveForm() {
  const e = litState.editing, row = { ...e.row };
  $("#litForm").querySelectorAll("[data-lit]").forEach((el) => { row[el.dataset.lit] = el.type === "checkbox" ? (el.checked ? "yes" : "no") : el.value; });
  if (!["name", "forward_primer", "reverse_primer", "probe", "citation"].some((f) => String(row[f] || "").trim())) {
    litState.errors = [{ field: "name", message: "입력된 내용이 없습니다 — 이름, 서열, 출처를 입력하세요" }]; e.row = row; litRenderForm(); return;
  }
  const rows = litState.payload.viruses[litState.virus].rows.map((r) => ({ ...r }));
  if (e.index === null) rows.push(row); else rows[e.index] = row;
  if (await litPostRows(rows)) { litState.editing = null; litState.errors = []; renderLiteraturePage(); $("#litStatus").textContent = "저장했습니다"; }
  else { e.row = row; litRenderForm(); }
}

async function litDelete(i) {
  const rows = litState.payload.viruses[litState.virus].rows;
  if (!await uiConfirm(`'${rows[i].name}' 항목을 목록에서 삭제할까요? (분석 목록에 이미 추가한 세트는 그대로 남습니다)`)) return;
  if (await litPostRows(rows.filter((_, k) => k !== i))) renderLiteraturePage();
}

async function litAddToAnalysis(i) {
  const r = litState.payload.viruses[litState.virus].rows[i];
  const warn = r.verified === "yes" ? "" : "\n\n⚠ 이 서열은 아직 원문과 대조하지 않았습니다. 결과는 서열이 맞다는 전제에서만 의미가 있습니다.";
  if (!await uiConfirm(`'${r.name}'을(를) 분석 목록(프라이머 입력 / 관리)에 추가할까요?${warn}`)) return;
  const res = await fetch("/api/literature/add", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ virus: litState.virus, lit_id: r.lit_id }) });
  const data = await res.json();
  if (!res.ok) { $("#litStatus").textContent = `추가 실패: ${data.error || (data.errors || []).map((x) => x.message).join(" / ")}`; return; }
  litState.payload = await fetchJSON("/api/literature");
  litState.virus = litState.payload.edition;
  renderLiteraturePage();
  $("#litStatus").textContent = `${data.name}: ${data.note}`;
}

$("#litBtn").addEventListener("click", () => showLiteraturePage());


// ---------------- 엔진 검증: Go Probe 설계 엔진 vs 파이썬(정답지) ----------------
const ENG = { sets: [], results: {} };
async function showEngineCheckPage() {
  state.activeId = null; state.activeType = null;
  renderSidebar();
  showMainView("engineView");
  const view = $("#engineView");
  view.innerHTML = "<p class='note'>불러오는 중...</p>";
  let st;
  try { st = await fetchJSON("/api/engine-check"); } catch (e) { view.innerHTML = `<p class="note">불러오기 실패: ${escHtml(e.message)}</p>`; return; }
  ENG.sets = st.sets || [];
  const selftest = (st.selftest || []).map((l) => `<div class="mono ${l.startsWith("PASS") ? "" : "bad"}" style="font-size:12px;">${escHtml(l)}</div>`).join("") || "<span class='note'>없음</span>";
  view.innerHTML = `
    <h2>🧪 엔진 검증 <span class="note">Probe 설계: Go 엔진 ↔ 파이썬(정답지)</span></h2>
    <div class="note" style="margin-bottom:10px;">같은 정렬(MSA)과 같은 조건을 두 엔진에 넣고, 후보 서열·가닥·위치·IUPAC 수·Tm/GC·커버리지·순서가 한 글자도 다르지 않은지 비교합니다. 규칙 기반 후보(own)만 비교하며 Primer3와 MFEprimer 검증은 파이썬에 남아 있습니다. 파이썬 결과가 정답입니다.</div>
    <div class="card" style="margin-bottom:14px;"><h3>Go 엔진 상태</h3>
      <div>${st.engine_built ? "✅ 엔진 준비됨" : "⚠️ 엔진을 빌드하지 못했습니다 (Go 설치 필요)"} <span class="note">${escHtml(st.go || "")}</span></div>
      <div style="margin-top:6px;"><b>자체 단위 검사</b> <span class="note">(파이썬이 계산한 값과 비트 단위로 비교)</span>${selftest}</div></div>
    <div class="card" style="margin-bottom:14px;"><h3>실제 세트로 비교</h3>
      <div style="display:flex;gap:10px;flex-wrap:wrap;align-items:end;">
        <label class="note">세트 <select id="engSet">${ENG.sets.map((s) => `<option value="${escHtml(s.set_id)}">${escHtml(s.name)} (${escHtml(s.set_id)})</option>`).join("")}</select></label>
        <label class="note">종류 <select id="engKind"><option value="design">Probe 설계</option><option value="scan">전체 유전체 스캔 (F/R·Probe 검색)</option></select></label>
        <label class="note">후보 수 <input id="engTop" type="number" min="1" max="60" value="10" style="width:60px;"></label>
        <label class="note">IUPAC 최대 <input id="engIupac" type="number" min="0" max="4" value="2" style="width:50px;"></label>
        <button id="engRun" class="btn btn-primary">선택한 세트 비교</button>
        <button id="engAll" class="btn btn-secondary" title="모든 세트를 IUPAC 0·2·4 조건으로 차례로 비교합니다">모든 세트 일괄 비교</button>
      </div>
      <div id="engOut" style="margin-top:12px;"></div></div>`;
  const out = view.querySelector("#engOut");
  const one = async (id, top, iup) => fetchJSON(`/api/engine-check/${encodeURIComponent(id)}?top=${top}&max_iupac=${iup}${view.querySelector("#engKind").value === "scan" ? "&kind=scan" : ""}`);
  view.querySelector("#engRun").addEventListener("click", async (e) => {
    const id = view.querySelector("#engSet").value;
    e.currentTarget.disabled = true; out.innerHTML = "<p class='note'>두 엔진을 실행하는 중...</p>";
    try { const r = await one(id, view.querySelector("#engTop").value, view.querySelector("#engIupac").value); out.innerHTML = engineResultHtml(r); }
    catch (err) { out.innerHTML = `<p class="note bad">실패: ${escHtml(err.message)}</p>`; }
    e.currentTarget.disabled = false;
  });
  view.querySelector("#engAll").addEventListener("click", async (e) => {
    const btn = e.currentTarget; btn.disabled = true;
    const rows = []; let bad = 0, py = 0, go = 0;
    const draw = (running) => { out.innerHTML = `<table class="data-table"><thead><tr><th>세트</th><th>IUPAC</th><th>후보</th><th>결과</th><th>Python 초</th><th>Go 초</th></tr></thead><tbody>${rows.join("")}</tbody></table>
      <p class="note" style="margin-top:8px;">${running ? "⏳ 진행 중... " : ""}${rows.length}건 · 불일치 ${bad}건 · 합계 Python ${py.toFixed(1)}초 / Go ${go.toFixed(1)}초</p>`; };
    const tid = "engine-all";
    let n = 0; const total = ENG.sets.length * (view.querySelector("#engKind").value === "scan" ? 1 : 3);
    const scanKind = view.querySelector("#engKind").value === "scan";
    for (const s of ENG.sets) for (const [iup, top] of (scanKind ? [[2, 10]] : [[0, 10], [2, 10], [4, 30]])) {
      toast(`엔진 일괄 비교 ${++n}/${total}`, { type: "run", id: tid });
      try {
        const r = await one(s.set_id, top, iup);
        const ok = r.ok; if (!ok) bad++; py += r.python_sec; go += r.go_sec;
        rows.push(`<tr><td>${escHtml(s.set_id)}</td><td>${iup}</td><td>${r.n_candidates}</td><td>${ok ? "✅ 동일" : "❌ " + escHtml(r.checks.filter((c) => !c.ok).map((c) => c.name).slice(0, 3).join(", "))}</td><td>${r.python_sec}</td><td>${r.go_sec}</td></tr>`);
      } catch (err) { bad++; rows.push(`<tr><td>${escHtml(s.set_id)}</td><td>${iup}</td><td></td><td>❌ ${escHtml(err.message)}</td><td></td><td></td></tr>`); }
      draw(true);
    }
    draw(false);
    toast(bad ? `엔진 비교: 불일치 ${bad}건` : `엔진 비교: ${rows.length}건 모두 동일`, { type: bad ? "err" : "ok", id: tid });
    btn.disabled = false;
  });
}
function engineResultHtml(r) {
  const head = r.ok ? `<div class="assessment-verdict ok">✅ 두 엔진의 결과가 완전히 같습니다</div>` : `<div class="assessment-verdict bad">❌ 다른 항목이 있습니다</div>`;
  const sp = r.go_sec ? (r.python_sec / r.go_sec) : 0;
  const rows = r.checks.map((c) => `<tr><td>${escHtml(c.name)}</td><td>${c.ok ? "✅" : "❌"}</td><td class="mono" style="font-size:12px;">${c.ok ? (typeof c.python === "string" ? escHtml(c.python) : "") : escHtml(JSON.stringify(c.python)) + "<br>↔ " + escHtml(JSON.stringify(c.go))}</td></tr>`).join("");
  return `${head}<div class="note" style="margin:6px 0;">유전체 ${r.n_rows}개 · 검사한 후보 ${r.n_candidates}개 · Python <b>${r.python_sec}초</b> · Go <b>${r.go_sec}초</b>${sp ? ` (Go가 ${sp >= 1 ? sp.toFixed(1) + "배 빠름" : (1 / sp).toFixed(1) + "배 느림"})` : ""}</div>
    <table class="data-table"><thead><tr><th>항목</th><th></th><th>값 / 차이</th></tr></thead><tbody>${rows}</tbody></table>`;
}
$("#engineBtn").addEventListener("click", showEngineCheckPage);
