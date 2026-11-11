// node --test tests/js/msa_column_map.test.mjs
// The pi bars above the MSA must sit on the same reference positions as the alignment columns below them.
import fs from "node:fs";
import path from "node:path";
import test from "node:test";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const src = fs.readFileSync(path.join(here, "../../webapp/static/app.js"), "utf8");
const body = src.slice(src.indexOf("function msaColumnToRefPos"), src.indexOf("// The pi bars are drawn in an SVG"));
const { msaColumnToRefPos, buildMsaRuler } = new Function(body + "; return { msaColumnToRefPos, buildMsaRuler };")();

const bind = (role, strand, start, end) => ({ role, strand, start: String(start), end: String(end), found: true });

test("F on the minus strand: columns run DOWN the reference coordinates", () => {
  const d = { msa: { columns: 158, consensus: "A".repeat(158), f_col: [0, 20], p_col: [54, 77], r_col: [138, 158] },
    binding: [bind("F", "-", 2482, 2501), bind("P", "-", 2420, 2442), bind("R", "+", 2344, 2363)] };
  const f = msaColumnToRefPos(d, 2344);
  assert.equal(f(0), 2501); assert.equal(f(19), 2482); assert.equal(f(138), 2363); assert.equal(f(157), 2344);
  assert.equal(f(54), 2442);
  for (let c = 1; c < 158; c++) assert.ok(f(c) < f(c - 1), `column ${c} must lie below column ${c - 1}`);
});

test("F on the plus strand: columns run UP (and the old start+column formula agrees)", () => {
  const d = { msa: { columns: 110, consensus: "A".repeat(110), f_col: [0, 20], p_col: [33, 64], r_col: [87, 110] },
    binding: [bind("F", "+", 117, 136), bind("P", "+", 150, 180), bind("R", "-", 204, 226)] };
  const f = msaColumnToRefPos(d, 117);
  for (const c of [0, 19, 33, 63, 87, 109]) assert.equal(f(c), 117 + c);
});

test("insertion columns (consensus gap) have no reference position and do not shift the rest", () => {
  const cons = "AA" + "A".repeat(20) + "-".repeat(5) + "A".repeat(20) + "AAA";       // 5 insertion columns after F
  const d = { msa: { columns: cons.length, consensus: cons, f_col: [2, 22], p_col: [27, 47], r_col: [47, 50] },
    binding: [bind("F", "-", 1359, 1378), bind("R", "+", 1336, 1338)] };      // 21 reference bases between F and R
  const f = msaColumnToRefPos(d, 1336);
  assert.equal(f(2), 1378); assert.equal(f(21), 1359);
  for (let c = 22; c < 27; c++) assert.equal(f(c), null);
  assert.equal(f(27), 1358);                                                           // the very next reference base
  assert.equal(f(49), 1336);
});

test("falls back to start+column when there is nothing to anchor on", () => {
  const f = msaColumnToRefPos({ msa: { columns: 10, consensus: "A".repeat(10) }, binding: [] }, 100);
  assert.equal(f(0), 100); assert.equal(f(9), 109);
  const g = msaColumnToRefPos({ msa: { columns: 10, consensus: "A".repeat(10), f_col: [0, 5] }, binding: [bind("F", "+", 100, 120)] }, 100);   // span disagrees
  assert.equal(g(3), 103);
});

test("the pi SVG is sized from the measured sequence width, not from a bare 'N x ch' assumption", () => {
  const upd = src.slice(src.indexOf("function updateMsaDiversityBaseline"), src.indexOf("function renderSetDiversitySection"));
  assert.match(upd, /fitMsaDiversitySvg\(\)/);
  assert.match(src, /function fitMsaDiversitySvg[\s\S]*getBoundingClientRect\(\)\.width[\s\S]*svg\.style\.width/);
});

test("ruler: one character per column, reference coordinates every 10, reading in the alignment's direction", () => {
  const d = { msa: { columns: 40, consensus: "A".repeat(40), f_col: [0, 20], r_col: [30, 40] },
    binding: [bind("F", "-", 2482, 2501), bind("R", "+", 2462, 2471)] };
  const ruler = buildMsaRuler(40, msaColumnToRefPos(d, 2462));
  assert.equal(ruler.length, 40);
  assert.equal(ruler.indexOf("|2500"), 1);                              // column 1 <-> position 2500 (F's 5' end is 2501 at column 0)
  assert.ok(ruler.indexOf("|2490") === 11 && ruler.indexOf("|2480") === 21 && ruler.indexOf("|2470") === 31);
  assert.ok(ruler.indexOf("|2500") < ruler.indexOf("|2490"));           // decreasing to the right
});

test("ruler: no tick starts on an insertion column, and the text keeps its length", () => {
  const cons = "AA" + "A".repeat(20) + "-".repeat(5) + "A".repeat(20) + "AAA";
  const d = { msa: { columns: cons.length, consensus: cons, f_col: [2, 22], r_col: [47, 50] },
    binding: [bind("F", "-", 1359, 1378), bind("R", "+", 1336, 1338)] };
  const ruler = buildMsaRuler(cons.length, msaColumnToRefPos(d, 1336));
  assert.equal(ruler.length, cons.length);
  assert.ok(!ruler.slice(22, 27).includes("|"));                        // no tick starts on an insertion column (a label may run over them)
  assert.equal(ruler.indexOf("|1360"), 20);                             // 1360 sits 2 columns before the insertion block
});

test("ruler without a position map is a plain tick every 10 columns", () => {
  assert.equal(buildMsaRuler(25, null), "|         |         |    ");
});

test("column selection measures text, not the padded box", () => {
  assert.match(src, /function msaCharWidth[\s\S]*createRange\(\)[\s\S]*function msaTextLeft/);
});
