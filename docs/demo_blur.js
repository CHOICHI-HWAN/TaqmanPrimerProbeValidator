// Demo mode for recording screenshots / GIFs of the dashboard: blurs everything that must not be shown.
//   - vendor / supplier names (see VENDOR_RE; edit the list for your own partners)
//   - nucleotide strings (primers, probes, standard materials, genome pieces): 10 or more bases in a row
//   - the primer/probe row, the consensus row and the amplicon alignment of the MSA / binding pop-ups
//   - sequence input boxes (primer forms, standard-material forms) and standard-material names in pop-ups
// Paste this file into the browser console (or run it with the browser tool) on the dashboard page; it keeps working while the
// page re-renders. `window.__demoBlurOff()` removes it.
(function () {
  const VENDOR_RE = /(?<![A-Za-z])(VendorA|VendorA|VendorB|VendorC|VendorD)(?![A-Za-z])/gi;
  const SEQ_RE = /(?<![A-Za-z0-9_])[ACGTURYKMSWBDHVNacgturykmswbdhvn]{10,}(?![A-Za-z0-9_])/g;
  const BLUR_SELECTORS = [".msa-arrow-row .msa-seq", ".msa-consensus .msa-seq", ".amp-align", ".standard-material-block .amp-align",
    "textarea", ".seq-box", ".mono.seq", "td.seq", ".primer-seq",
    ".standard-material-id", ".modal-box table td:first-child"].join(",");        // standard-material names (pop-ups)
  const PX = 7;

  // extra words to hide anywhere they appear (e.g. the names of your standard materials): set window.__demoBlurWords = [...] BEFORE
  // loading this file; the dashboard's standard-material names are added automatically below.
  // pure: [start, end) ranges of the text that must be hidden (exported for the tests)
  function sensitiveRanges(text, words) {
    const out = [];
    const res = [VENDOR_RE, SEQ_RE];
    if (words && words.length) res.push(new RegExp("(?<![A-Za-z0-9])(" + words.map((x) => x.replace(/[.*+?^${}()|[\]\\]/g, "\\$&")).join("|") + ")(?![A-Za-z0-9])", "g"));
    for (const re of res) {
      re.lastIndex = 0;
      let m;
      while ((m = re.exec(text))) {
        const run = m[0];
        if (re === SEQ_RE && (new Set(run.toUpperCase()).size === run.length || new Set(run.toUpperCase()).size === 1)) continue;   // alphabets / NNNN
        out.push([m.index, m.index + run.length]);
      }
    }
    out.sort((a, b) => a[0] - b[0]);
    const merged = [];
    for (const r of out) {
      const last = merged[merged.length - 1];
      if (last && r[0] <= last[1]) last[1] = Math.max(last[1], r[1]); else merged.push(r.slice());
    }
    return merged;
  }
  if (typeof module !== "undefined") { module.exports = { sensitiveRanges }; }
  if (typeof document === "undefined") return;

  const style = document.createElement("style");
  style.id = "demoBlurStyle";
  style.textContent = `.demo-blur { filter: blur(${PX}px); user-select: none; }
    textarea, input[type=text]:not(.mt-menu-search) { filter: blur(${PX}px); }`;
  document.head.appendChild(style);

  function wrapText(node) {
    const text = node.nodeValue;
    const ranges = sensitiveRanges(text, window.__demoBlurWords);
    if (!ranges.length) return;
    const frag = document.createDocumentFragment();
    let pos = 0;
    for (const [a, b] of ranges) {
      if (a > pos) frag.appendChild(document.createTextNode(text.slice(pos, a)));
      const s = document.createElement("span");
      s.className = "demo-blur";
      s.textContent = text.slice(a, b);
      frag.appendChild(s);
      pos = b;
    }
    if (pos < text.length) frag.appendChild(document.createTextNode(text.slice(pos)));
    node.parentNode.replaceChild(frag, node);
  }

  function sweep(root) {
    const walker = document.createTreeWalker(root, NodeFilter.SHOW_TEXT, {
      acceptNode: (n) => (n.parentElement && !n.parentElement.closest("script,style,.demo-blur") && n.nodeValue.length >= 3) ? NodeFilter.FILTER_ACCEPT : NodeFilter.FILTER_REJECT,
    });
    const nodes = [];
    while (walker.nextNode()) nodes.push(walker.currentNode);
    nodes.forEach(wrapText);
    root.querySelectorAll(BLUR_SELECTORS).forEach((e) => { e.classList.add("demo-blur"); });
    // tooltips / titles / placeholders can carry names too
    root.querySelectorAll("[title],[placeholder]").forEach((e) => {
      for (const a of ["title", "placeholder"]) if (e.hasAttribute(a) && sensitiveRanges(e.getAttribute(a), window.__demoBlurWords).length) e.removeAttribute(a);
    });
  }

  let timer = null;
  const obs = new MutationObserver(() => { clearTimeout(timer); timer = setTimeout(() => { obs.disconnect(); sweep(document.body); obs.observe(document.body, { childList: true, subtree: true, characterData: true }); }, 60); });
  sweep(document.body);
  obs.observe(document.body, { childList: true, subtree: true, characterData: true });
  window.__demoBlurOff = () => { obs.disconnect(); style.remove(); document.querySelectorAll(".demo-blur").forEach((e) => e.classList.remove("demo-blur")); };
})();
