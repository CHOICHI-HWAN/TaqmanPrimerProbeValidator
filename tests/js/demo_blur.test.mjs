// node --test tests/js/demo_blur.test.mjs -- the recording helper hides vendor names and nucleotide strings, nothing else.
import { createRequire } from "node:module";
import path from "node:path";
import test from "node:test";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";
import { synth } from "./seqjs.mjs";

const here = path.dirname(fileURLToPath(import.meta.url));
const { sensitiveRanges } = createRequire(import.meta.url)(path.join(here, "../../docs/demo_blur.js"));
const cut = (t, words) => sensitiveRanges(t, words).map(([a, b]) => t.slice(a, b));

test("listed words (vendor names) are hidden wherever they sit in a name, never inside longer words", () => {
  const w = ["VendorA", "VendorB"];
  assert.deepEqual(cut("LASV-VendorA-S-GPC-PREV1-1B", w), ["VendorA"]);
  assert.deepEqual(cut("VendorA_L_1 and VendorB_L_1 / VendorB", w), ["VendorA", "VendorB", "VendorB"]);
  assert.deepEqual(cut("VendorAB, xVendorA, VendorAx", w), []);
  assert.deepEqual(cut("LASV-VendorA-S", []), []);                       // nothing is built in
});

test("nucleotide strings of 10+ bases are hidden, short ones and alphabets are not", () => {
  const s = synth(24, 3);
  assert.deepEqual(cut(`primer ${s} ok`), [s]);
  assert.deepEqual(cut("ACGT RYKM 123456"), []);
  assert.deepEqual(cut("ABCDGHKMNRSTVWY"), []);                          // an IUPAC alphabet, every letter different
  assert.deepEqual(cut("NNNNNNNNNNNN"), []);
});

test("a name next to a sequence gives two separate ranges, in order", () => {
  const s = synth(12, 5);
  const text = `VendorA-${s}`;
  assert.deepEqual(sensitiveRanges(text, ["VendorA"]), [[0, 7], [8, 20]]);
});

test("extra words (standard-material names) are hidden as whole words only", () => {
  const cutW = (t, w) => sensitiveRanges(t, w).map(([a, b]) => t.slice(a, b));
  assert.deepEqual(cutW("매치 O (GC2) and GC22, xGC2", ["GC2", "NP1"]), ["GC2"]);
  assert.deepEqual(cutW("none here", ["GC2"]), []);
  assert.deepEqual(cutW("a.b (x) [y] aXb", ["a.b", "(x)"]), ["a.b", "(x)"]);   // regex characters in a name are taken literally
});
