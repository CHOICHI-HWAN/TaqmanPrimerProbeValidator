// "Firebase 업로드" button for the local analysis dashboard (top right). Crawls this dashboard's own read-only GET endpoints,
// scrubs server paths, and writes the result to Firestore (results/<key>~<chunk>) as the signed-in admin. Analysis stays local.
// Pure helpers are exported for the node tests; everything that touches the page is behind `typeof document`.
const SDK = "https://www.gstatic.com/firebasejs/10.12.2";
export const CHUNK = 600 * 1024;
// Endpoints whose content depends on the host scope are exported once per scope (key suffix "@human" etc.; no suffix = the server's
// default scope). The read-only site picks the variant from the ?host= the dashboard sends. Everything else is exported once.
export const SCOPES = ["human", "animal", "unknown", "all"];
export const SCOPED_RE = /^(?:(?:sets|groups)(?:\/[^/]+)?|overview|host-scope)$/;
const FIXED = ["profile", "diversity", "metadata", "metadata/table", "standards/representativeness",
  "standards/manifest", "primers/manifest", "changelog", "requirements", "params", "jobs/active", "pipeline/freshness", "ui-state"];
const FIXED_SCOPED = ["overview", "sets", "groups", "host-scope"];
const IMG_CT = new Set(["image/png", "image/jpeg", "image/gif", "image/webp"]);

export function scrub(v) {                                   // absolute server paths -> "…/<last part>"
  if (typeof v === "string") return /^\/(home|BI|mnt|tmp)\//.test(v) ? "…/" + v.replace(/\/+$/, "").split("/").pop() : v;
  if (Array.isArray(v)) return v.map(scrub);
  if (v && typeof v === "object") return Object.fromEntries(Object.entries(v).map(([k, x]) => [k, scrub(x)]));
  return v;
}

export const docId = (key) => key.replaceAll("/", "~");

export function b64(bytes) {
  let s = "";
  for (let i = 0; i < bytes.length; i += 0x8000) s += String.fromCharCode.apply(null, bytes.subarray(i, i + 0x8000));
  return btoa(s);
}

export function chunks(bytes) {                              // -> [Uint8Array] of at most CHUNK bytes (at least one)
  const n = Math.max(1, Math.ceil(bytes.length / CHUNK));
  return Array.from({ length: n }, (_, i) => bytes.subarray(i * CHUNK, (i + 1) * CHUNK));
}

async function gzip(bytes) {
  return new Uint8Array(await new Response(new Blob([bytes]).stream().pipeThrough(new CompressionStream("gzip"))).arrayBuffer());
}

/** Everything to publish: [{key, gz, ctype, bytes}] (bytes already gzip-compressed when gz). `get(path)` -> Response. */
export async function collect(project, get, progress = () => {}) {
  const out = [], enc = new TextEncoder();
  const addJson = async (key, obj, ctype = "application/json; charset=utf-8") => {
    out.push({ key: `${project}/${key}`, gz: true, ctype, bytes: await gzip(enc.encode(typeof obj === "string" ? obj : JSON.stringify(scrub(obj)))) });
  };
  const json = async (path, scope) => {
    const r = await get("/api/" + path + (scope ? (path.includes("?") ? "&" : "?") + "host=" + scope : ""));
    return r.ok ? r.json() : null;
  };
  const data = {};
  for (const p of FIXED) {
    progress(`읽는 중: ${p}`);
    data[p] = await json(p);
    if (data[p] !== null) await addJson(p, data[p]);
  }
  for (const scope of [null, ...SCOPES]) {                          // null = default scope, then each scope explicitly
    const sfx = scope ? `@${scope}` : "";
    for (const p of FIXED_SCOPED) {
      progress(`읽는 중: ${p}${sfx}`);
      const d = await json(p, scope);
      if (d !== null) { await addJson(p + sfx, d); if (!scope) data[p] = d; }
    }
    for (const s of data.sets || []) { progress(`세트 ${s.set_id}${sfx}`); const d = await json(`sets/${encodeURIComponent(s.set_id)}`, scope); if (d) await addJson(`sets/${s.set_id}${sfx}`, d); }
    for (const g of data.groups || []) { progress(`그룹 ${g.group_id}${sfx}`); const d = await json(`groups/${encodeURIComponent(g.group_id)}`, scope); if (d) await addJson(`groups/${g.group_id}${sfx}`, d); }
  }
  for (const seg of Object.keys((data.profile || {}).segments || {})) {
    const u = await get(`/api/standards/unmatched?seg=${seg}`);
    if (u.ok) out.push({ key: `${project}/standards/unmatched/${seg}`, gz: true, ctype: "text/plain; charset=utf-8", bytes: await gzip(new Uint8Array(await u.arrayBuffer())) });
    const a = await get(`/api/auspice/${seg}`);
    if (a.ok) { try { await addJson(`auspice/${seg}`, await a.json()); } catch { /* not a tree JSON */ } }
  }
  const ui = data["ui-state"] || {};
  const image = async (path, key) => {
    const r = await get(path);
    const ct = (r.headers.get("Content-Type") || "").split(";")[0].trim();
    if (r.ok && IMG_CT.has(ct)) out.push({ key: `${project}/${key}`, gz: false, ctype: ct, bytes: new Uint8Array(await r.arrayBuffer()) });
  };
  for (const [ent, imgs] of Object.entries(ui.note_images || {})) {
    const kind = ent.startsWith("G") ? "groups" : "sets";
    for (const im of imgs) await image(`/api/${kind}/${ent}/note/images/${im.id}`, `${kind}/${ent}/note/images/${im.id}`);
  }
  for (const [ent, vers] of Object.entries(ui.verifications || {})) {
    const kind = ent.startsWith("G") ? "groups" : "sets";
    for (const v of vers) for (const im of v.images || []) await image(`/api/${kind}/${ent}/verify/${v.id}/images/${im.id}`, `${kind}/${ent}/verify/${v.id}/images/${im.id}`);
  }
  return out;
}

export function explainError(e) {
  const code = (e && e.code) || "";
  const msg = (e && e.message) || String(e);
  if (code === "permission-denied") return "관리자 계정(admin_emails.txt)으로 로그인해야 올릴 수 있습니다.";
  if (code === "auth/popup-blocked") return "브라우저가 로그인 팝업을 막았습니다. 주소창의 팝업 차단 표시에서 이 사이트(localhost)의 팝업을 허용한 뒤 다시 누르세요.";
  if (code === "auth/popup-closed-by-user" || code === "auth/cancelled-popup-request") return "로그인 창이 닫혀서 취소되었습니다.";
  if (code === "auth/unauthorized-domain") return "Firebase 콘솔 > Authentication > 설정 > 승인된 도메인에 이 주소의 호스트(localhost)를 추가해야 합니다.";
  if (["auth/network-request-failed", "auth/internal-error", "auth/web-storage-unsupported", "auth/timeout"].includes(code))
    return "로그인 연동이 막혔습니다. Firefox는 주소창 왼쪽 방패 아이콘에서 이 사이트의 '향상된 추적 방지'를 끄면 되는 경우가 많습니다(그래도 안 되면 Chrome 사용). (" + code + ")";
  return code ? `${msg} (${code})` : msg;
}

if (typeof document !== "undefined") {
  const btn = document.createElement("button");
  btn.id = "fbPublishBtn";
  btn.textContent = "☁ Firebase 업로드";
  btn.hidden = true;
  btn.style.cssText = "position:fixed;top:10px;right:16px;z-index:5000;padding:7px 14px;border-radius:8px;border:1px solid #94a3b8;background:#fff;color:#0f172a;cursor:pointer;font-size:13px;box-shadow:0 1px 4px rgba(0,0,0,.2)";
  const msg = document.createElement("div");                       // errors stay on screen (alert() alone is easy to miss)
  msg.id = "fbPublishMsg";
  msg.hidden = true;
  msg.style.cssText = "position:fixed;top:46px;right:16px;z-index:5000;max-width:380px;padding:8px 12px;border-radius:8px;font-size:12.5px;line-height:1.45;box-shadow:0 1px 6px rgba(0,0,0,.25);background:#fff;color:#0f172a;border:1px solid #94a3b8;cursor:pointer";
  msg.title = "클릭하면 닫힙니다";
  msg.onclick = () => { msg.hidden = true; };
  document.body.appendChild(btn);
  document.body.appendChild(msg);
  // step log: shows how far the button got (config -> libraries -> auth warm-up -> click), so a failure can be read off the screen
  const steps = document.createElement("div");
  steps.id = "fbPublishSteps";
  steps.style.cssText = "position:fixed;top:44px;right:16px;z-index:4999;font:11px/1.35 ui-monospace,monospace;color:#64748b;background:rgba(255,255,255,.85);padding:2px 6px;border-radius:6px;max-width:380px;text-align:right";
  steps.hidden = true;
  document.body.appendChild(steps);
  const stepLog = [];
  const step = (t) => { stepLog.push(t); steps.textContent = stepLog.slice(-4).join(" › "); steps.hidden = false; msg.style.top = "66px"; };
  window.addEventListener("error", (ev) => { if (String(ev.filename || "").includes("publish.js")) step("JS 오류: " + ev.message); });
  window.addEventListener("unhandledrejection", (ev) => { const m = String((ev.reason && (ev.reason.code || ev.reason.message)) || ev.reason); step("오류: " + m); });
  btn.addEventListener("click", () => step("클릭됨"), true);          // registered first: proves the click reaches the button
  step("스크립트 로드됨");

  const setLabel = (t) => { btn.textContent = t; };
  const say = (t, bad) => { msg.textContent = t; msg.style.borderColor = bad ? "#e5484d" : "#94a3b8"; msg.hidden = !t; };
  let sdk = null, busy = false;

  // The Firebase libraries are loaded and initialised when the page opens, NOT on click: the login popup must be opened inside the
  // click itself, or the browser (Firefox in particular) treats it as unsolicited and blocks it.
  async function preload(cfg) {
    const [{ initializeApp }, { getAuth, signInWithPopup, GoogleAuthProvider, getRedirectResult }, { getFirestore, doc, setDoc }] = await Promise.all(
      ["app", "auth", "firestore"].map((m) => import(`${SDK}/firebase-${m}.js`)));
    const app = initializeApp(cfg.firebaseConfig);
    sdk = { auth: getAuth(app), db: getFirestore(app), signInWithPopup, GoogleAuthProvider, doc, setDoc };
    // Warm-up: the first sign-in otherwise opens its popup ~0.5 s AFTER the click (it first loads Google's auth iframe), which is
    // past Firefox's ~1 s window for click-initiated popups on a slow connection. getRedirectResult() initialises that machinery now.
    getRedirectResult(sdk.auth).catch(() => null);
  }

  async function upload(cfg) {
    const items = await collect(cfg.project, (p) => fetch(p), (t) => setLabel(t));
    items.push({ key: "_probe", gz: false, ctype: "application/json", bytes: new TextEncoder().encode('{"ok":1}') });
    let done = 0;
    for (const it of items) {
      const parts = chunks(it.bytes);
      for (let i = 0; i < parts.length; i++) {
        const payload = { n: parts.length, gz: it.gz, ctype: it.ctype, data: b64(parts[i]) };
        for (let t = 1; ; t++) {
          try { await sdk.setDoc(sdk.doc(sdk.db, "results", `${docId(it.key)}~${i}`), payload); break; }
          catch (e) { if (t >= 3) throw e; await new Promise((r) => setTimeout(r, 1000 * t)); }
        }
      }
      setLabel(`업로드 ${++done}/${items.length}`);
    }
    return items.length;
  }

  fetch("/api/firebase/config").then((r) => { step("설정 응답 " + r.status); return r.ok ? r.json() : null; }).then((cfg) => {
    if (!cfg || !cfg.firebaseConfig || !cfg.project) { step("설정 없음(config/firebase.json)"); return; }
    btn.hidden = false;
    step("설정 읽음");
    setLabel("☁ Firebase 업로드 (준비 중…)");
    preload(cfg).then(() => { setLabel("☁ Firebase 업로드"); step("라이브러리·인증 준비됨"); }, (e) => {
      step("라이브러리 로드 실패: " + (e.message || e));
      console.error("Firebase SDK load failed", e);
      setLabel("☁ Firebase 업로드");
      say("Firebase 라이브러리를 불러오지 못했습니다(인터넷 연결이나 추적 차단 확장프로그램을 확인하세요): " + (e.message || e), true);
    });

    btn.onclick = () => {
      if (busy) return;
      if (!sdk) { say("Firebase 라이브러리를 아직 불러오는 중입니다. 잠시 후 다시 누르세요.", false); return; }
      say("", false);
      step(sdk.auth.currentUser ? "로그인됨" : "로그인 팝업 요청");
      busy = true; btn.disabled = true;
      // Called directly inside the click -- no await / confirm() before it.
      const login = sdk.auth.currentUser ? Promise.resolve() : (setLabel("로그인 중…"), sdk.signInWithPopup(sdk.auth, new sdk.GoogleAuthProvider()));
      // Firefox (tracking protection / partitioned storage) can let the popup log in and close without the result ever reaching this
      // page, so the promise just never settles. Say so after a while instead of leaving the user with a silent button.
      const stuck = sdk.auth.currentUser ? null : setTimeout(() => say("로그인 창에서 로그인했는데 여기서 진행이 안 되면, 브라우저가 로그인 결과 전달을 막은 것입니다. Firefox: 주소창 왼쪽 방패 아이콘 > '향상된 추적 방지'를 이 사이트에서 끄고 페이지를 새로고침한 뒤 다시 누르세요. 그래도 안 되면 Chrome에서 사용하세요.", true), 20000);
      login.then(async () => {
        if (stuck) clearTimeout(stuck);
        say("", false);
        if (!confirm(`${sdk.auth.currentUser.email} 계정으로 현재 분석 결과(분석 제외, 결과만)를 Firebase에 올립니다.\n허가된 계정만 열람할 수 있습니다. 계속할까요?`)) { setLabel("☁ Firebase 업로드"); return; }
        const n = await upload(cfg);
        setLabel("☁ Firebase 업로드");
        say(`업로드 완료: ${n}개 파일. https://primervalidator.web.app/${cfg.project}/ 에서 확인하세요.`, false);
      }).catch((e) => {
        if (stuck) clearTimeout(stuck);
        console.error("Firebase upload failed", e);
        setLabel("☁ Firebase 업로드");
        say("업로드 실패: " + explainError(e), true);
      }).finally(() => { busy = false; btn.disabled = false; });
    };
  }).catch((e) => { step("설정 요청 실패: " + (e.message || e)); });
}
