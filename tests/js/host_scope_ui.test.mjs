// node --test tests/js/host_scope_ui.test.mjs
// The host-scope toggle: every API read carries ?host=..., the statement above a result says what is counted.
import fs from "node:fs";
import path from "node:path";
import test from "node:test";
import assert from "node:assert/strict";
import { fileURLToPath } from "node:url";

const here = path.dirname(fileURLToPath(import.meta.url));
const src = fs.readFileSync(path.join(here, "../../webapp/static/app.js"), "utf8");
const block = src.slice(src.indexOf("// ---------------- host scope"), src.indexOf("function hostScopeSegments"));
const badgeSrc = src.slice(src.indexOf("function hostScopeBadgeHtml"), src.indexOf("\n}\n", src.indexOf("function hostScopeBadgeHtml")) + 3);

function setup(hostScope) {
  const calls = [];
  const window = { fetch: (url, init) => { calls.push(url); return Promise.resolve(); } };
  const state = {};
  const localStorage = { getItem: () => hostScope };
  new Function("state", "window", "localStorage", block)(state, window, localStorage);
  return { window, state, calls };
}

test("API reads carry the chosen host scope, other requests are left alone", async () => {
  const { window, state, calls } = setup("animal");
  assert.equal(state.hostScope, "animal");
  await window.fetch("/api/sets");
  await window.fetch("/api/sets/P001?x=1");
  await window.fetch("/api/sets?host=all");                       // an explicit host wins
  await window.fetch("/api/sets", { method: "POST" });            // writes are not touched
  await window.fetch("/vendor/d3.min.js");
  assert.deepEqual(calls, ["/api/sets?host=animal", "/api/sets/P001?x=1&host=animal", "/api/sets?host=all", "/api/sets", "/vendor/d3.min.js"]);
});

test("without a choice nothing is added (the server then uses its primary scope)", async () => {
  const { window, calls } = setup(null);
  await window.fetch("/api/overview");
  assert.deepEqual(calls, ["/api/overview"]);
});

test("the statement above a result names the scope and what is left out", () => {
  const HOST_SCOPE_LABEL = { human: "사람", animal: "동물", unknown: "미상", all: "전체" };
  const HOST_CLASS_LABEL = { human: "사람", animal: "동물", unknown: "미상", lab: "실험 숙주" };
  const f = new Function("HOST_SCOPE_LABEL", "HOST_CLASS_LABEL", "escHtml", badgeSrc + "; return hostScopeBadgeHtml;")(HOST_SCOPE_LABEL, HOST_CLASS_LABEL, (x) => x);
  const info = { S: { scope: "human", scope_counts: { human: 571, animal: 52, unknown: 32, all: 655 }, counts: { human: 571, animal: 52, unknown: 32, lab: 7 }, excluded_lab: 7 } };
  const html = f(info);
  assert.match(html, /숙주 기준: 사람/);
  assert.match(html, /571/);
  assert.match(html, /동물 52/);
  assert.match(html, /실험 숙주 7개/);
  assert.equal(f({}), "");
  info.S.scope = "all";
  assert.match(f(info), /사람·동물·미상 전체 게놈 S 655|전체 게놈 655/);
});
