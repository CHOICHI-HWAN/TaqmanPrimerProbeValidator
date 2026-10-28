// node --test tests/js/binding_status.test.mjs
// "Not validated" is not the same as "not bound": a genome whose primer is only over the tolerance, away from the 3' end, is a
// caution (yellow), not "미결합" (red). Probe mismatches and 3'-end mismatches stay red.
import fs from "node:fs";
import path from "node:path";
import test from "node:test";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const src = fs.readFileSync(path.join(here, "../../webapp/static/app.js"), "utf8");
const grab = (start) => { const i = src.indexOf(start); let d = 0; for (let k = src.indexOf("{", i); k < src.length; k++) { if (src[k] === "{") d++; if (src[k] === "}" && !--d) return src.slice(i, k + 1); } };
const code = ["const PASS = { mm: 1, pm: 0 };", grab("function verdictFor"), src.slice(src.indexOf("const PARTIAL_MM"), src.indexOf(";", src.indexOf("const PARTIAL_MM")) + 1),
  grab("function standardState"), grab("function nearBinding"), grab("function verdictClass")].join("\n");
const { verdictFor, nearBinding, verdictClass, standardState } = new Function(code + "; return { verdictFor, nearBinding, verdictClass, standardState };")();
const g = (f, p, r, extra = {}) => ({ f_mismatch: f, p_mismatch: p, r_mismatch: r, f_3p_ok: true, r_3p_ok: true, ...extra });

test("validated genome is ok", () => assert.equal(verdictClass(g(0, 0, 1), verdictFor(g(0, 0, 1))), "ok"));

test("primer over the pass tolerance but bound, probe perfect, 3' ends fine -> caution (yellow), not unbound", () => {
  const x = g(0, 0, 2);                                       // the Nipah MW535746 case
  assert.equal(verdictFor(x).ok, false);
  assert.equal(nearBinding(x), true);
  assert.equal(verdictClass(x, verdictFor(x)), "warn");
});

test("probe mismatch stays a failure", () => {
  const x = g(0, 1, 1);
  assert.equal(nearBinding(x), false);
  assert.equal(verdictClass(x, verdictFor(x)), "bad");
});

test("a mismatch at a primer's 3' end stays a failure", () => {
  assert.equal(nearBinding(g(0, 0, 1, { r_3p_ok: false })), false);
  assert.equal(nearBinding(g(2, 0, 0, { f_3p_ok: false })), false);
});

test("too many mismatches, or no binding at all, is a failure", () => {
  assert.equal(nearBinding(g(0, 0, 3)), false);
  assert.equal(nearBinding(g(null, 8, null)), false);
  assert.equal(nearBinding(g(1, null, 1)), false);
});

test("standard-material partial-match class is unchanged", () => {
  assert.equal(standardState(g(0, 0, 1)), "ok");
  assert.equal(standardState(g(0, 0, 2)), "partial");
  assert.equal(standardState(g(null, 6, null)), "none");
});
