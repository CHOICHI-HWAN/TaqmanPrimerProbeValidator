// Data backup / restore: a small icon button at the top-right of the nav bar opens a modal dialog with a progress bar.
// Self-contained (no app.js changes). Server side: /api/backup/list|download (GET), /api/backup/create|restore (POST) -> scripts/backup_data.py.
(function () {
  "use strict";
  const fmt = (n) => (n >= 1073741824 ? (n / 1073741824).toFixed(1) + " GB" : n >= 1048576 ? (n / 1048576).toFixed(1) + " MB" : Math.max(1, Math.round(n / 1024)) + " KB");
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  let dlg = null;

  function ensureStyle() {
    if (document.getElementById("bkStyle")) return;
    const st = document.createElement("style");
    st.id = "bkStyle";
    st.textContent = `
      .bk-icon-btn { margin-left:auto; width:30px; height:30px; border-radius:50%; border:1px solid var(--border,#ccc); background:var(--bg,#fff); color:var(--fg,#222);
        font-size:15px; line-height:1; cursor:pointer; display:inline-flex; align-items:center; justify-content:center; flex:none; }
      .bk-icon-btn:hover { filter:brightness(1.08); }
      .bk-progress-track { height:8px; border-radius:5px; background:var(--border,#ddd); overflow:hidden; margin:4px 0 10px; }
      .bk-progress-fill { height:100%; width:0%; background:var(--accent,#4C78A8); transition:width .3s; }
      .bk-progress-fill.indeterminate { width:40% !important; animation: bkSlide 1.1s ease-in-out infinite; }
      .bk-progress-fill.done { background:var(--ok,#54A24B); width:100% !important; }
      .bk-progress-fill.error { background:var(--bad,#E45756); width:100% !important; }
      @keyframes bkSlide { 0% { margin-left:-40%; } 100% { margin-left:100%; } }`;
    document.head.appendChild(st);
  }

  function ensureDialog() {
    if (dlg) return dlg;
    ensureStyle();
    dlg = document.createElement("div");
    dlg.id = "backupDialog";
    dlg.hidden = true;
    dlg.style.cssText = "position:fixed;inset:0;background:rgba(0,0,0,.45);z-index:9000;display:flex;align-items:center;justify-content:center;";
    dlg.innerHTML = `
      <div style="background:var(--panel,#fff);color:var(--text,#222);border:1px solid var(--border,#ccc);border-radius:10px;max-width:760px;width:92%;max-height:88vh;overflow:auto;padding:18px;">
        <div style="display:flex;justify-content:space-between;align-items:center;margin-bottom:8px;">
          <h3 style="margin:0;">💾 데이터 백업 · 복원</h3><button class="btn btn-secondary" id="bkClose">닫기</button></div>
        <div class="note" style="margin-bottom:10px;">이 에디션의 설정·입력 프라이머·프로필·다운로드 데이터·분석 결과를 하나의 압축 파일(.tar.xz, 빠르게는 .tar.zst)로 저장합니다. 코드는 포함하지 않습니다(git에 있음).
          프라이머 서열과 업체명이 들어 있어 <b>대외비</b>입니다. 저장 위치는 git 밖입니다.</div>
        <div style="display:flex;gap:14px;align-items:center;flex-wrap:wrap;margin-bottom:10px;">
          <label><input type="checkbox" id="bkData" checked> data 포함 (원본 서열, 클수록 오래 걸림)</label>
          <label><input type="checkbox" id="bkOutput" checked> output 포함 (분석 결과)</label>
          <label title="zstd -10: 약 10배 빠르지만 파일이 약 1.7배 큽니다. 끄면 xz(가장 작음)"><input type="checkbox" id="bkFast"> 빠르게 (zstd, 용량 큼)</label>
          <button class="btn btn-primary" id="bkCreate">백업 만들기</button></div>
        <div id="bkStatus" class="note" style="min-height:1.4em;"></div>
        <div class="bk-progress-track" id="bkProgressTrack" hidden><div class="bk-progress-fill" id="bkProgressFill"></div></div>
        <table class="data-table" style="width:100%;"><thead><tr><th>만든 시각</th><th>버전</th><th>크기(원본)</th><th>포함</th><th></th></tr></thead><tbody id="bkRows"></tbody></table>
      </div>`;
    document.body.appendChild(dlg);
    dlg.addEventListener("click", (e) => { if (e.target === dlg) close(); });
    dlg.querySelector("#bkClose").onclick = close;
    dlg.querySelector("#bkCreate").onclick = create;
    return dlg;
  }

  function close() { if (dlg) dlg.hidden = true; }
  function status(msg, bad) { const el = dlg.querySelector("#bkStatus"); el.textContent = msg; el.style.color = bad ? "var(--bad,#c0392b)" : ""; }
  function progress(state) {
    // state: "running" | "done" | "error" | null (hide)
    const track = dlg.querySelector("#bkProgressTrack"), fill = dlg.querySelector("#bkProgressFill");
    if (!state) { track.hidden = true; return; }
    track.hidden = false;
    fill.className = "bk-progress-fill" + (state === "running" ? " indeterminate" : state === "done" ? " done" : state === "error" ? " error" : "");
  }

  async function refresh() {
    const r = await fetch("/api/backup/list").then((x) => x.json()).catch(() => null);
    const body = dlg.querySelector("#bkRows");
    if (!r) { body.innerHTML = "<tr><td colspan=5>목록을 읽지 못했습니다.</td></tr>"; return; }
    body.innerHTML = r.items.length ? r.items.map((it) => `<tr>
      <td>${esc(it.created.replace("T", " "))}</td><td>${esc(it.version)}</td><td>${fmt(it.archive_bytes)} <span class="note">(${fmt(it.uncompressed_bytes)})</span></td>
      <td class="note">${esc(it.includes.join(", "))}</td>
      <td style="white-space:nowrap;"><a class="btn btn-secondary" href="/api/backup/download?file=${encodeURIComponent(it.file)}" download>내려받기</a>
        <button class="btn btn-secondary" data-restore="${esc(it.file)}">복원</button></td></tr>`).join("")
      : "<tr><td colspan=5 class='note'>아직 백업이 없습니다.</td></tr>";
    body.querySelectorAll("[data-restore]").forEach((b) => (b.onclick = () => restore(b.dataset.restore)));
  }

  async function watch(jobId, doneText) {
    progress("running");
    for (;;) {
      await new Promise((r) => setTimeout(r, 1500));
      const j = await fetch("/api/jobs/" + jobId).then((x) => x.json()).catch(() => null);
      if (!j) continue;
      status(j.progress || "…", j.status === "error");
      if (j.status !== "running") {
        progress(j.status === "done" ? "done" : "error");
        if (j.status === "done") status(doneText + " " + (j.progress || ""));
        break;
      }
    }
    dlg.querySelectorAll("button").forEach((b) => (b.disabled = false));
    refresh();
  }

  async function post(url, payload) {
    const res = await fetch(url, { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(payload) });
    const j = await res.json().catch(() => ({}));
    if (!res.ok) throw new Error(j.error || res.statusText);
    return j;
  }

  async function create() {
    try {
      dlg.querySelectorAll("button").forEach((b) => (b.disabled = true));
      const j = await post("/api/backup/create", { data: dlg.querySelector("#bkData").checked, output: dlg.querySelector("#bkOutput").checked, fast: dlg.querySelector("#bkFast").checked });
      status("백업 중…");
      watch(j.job_id, "백업 끝.");
    } catch (e) { status(String(e.message || e), true); progress("error"); dlg.querySelectorAll("button").forEach((b) => (b.disabled = false)); }
  }

  async function restore(file) {
    if (!confirm(`${file}\n\n이 백업으로 설정·입력·프로필·data·output을 되돌립니다.\n지금 폴더는 지우지 않고 .pre_restore_* 폴더로 옮겨 둡니다. 계속할까요?`)) return;
    try {
      dlg.querySelectorAll("button").forEach((b) => (b.disabled = true));
      const j = await post("/api/backup/restore", { file, confirm: true });
      status("복원 중…");
      watch(j.job_id, "복원 끝. 화면을 새로고침하세요.");
    } catch (e) { status(String(e.message || e), true); progress("error"); dlg.querySelectorAll("button").forEach((b) => (b.disabled = false)); }
  }

  function open() { ensureDialog(); dlg.hidden = false; status(""); progress(null); refresh(); }

  function install() {
    const bar = document.getElementById("navBar");
    if (!bar || document.getElementById("backupIconBtn")) return;
    ensureStyle();
    const b = document.createElement("button");
    b.id = "backupIconBtn";
    b.className = "bk-icon-btn";
    b.type = "button";
    b.title = "💾 데이터 백업 · 복원";
    b.textContent = "💾";
    b.onclick = open;
    bar.appendChild(b);
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", install); else install();
})();
