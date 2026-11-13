// Synthetic nucleotide strings for JS tests (the repository holds no real sequences; see tests/seqgen.py).
export function synth(n, seed) {
  let x = seed * 2654435761 % 4294967296;
  const next = () => { x = (x * 1664525 + 1013904223) % 4294967296; return x / 4294967296; };
  const out = [];
  while (out.length < n) {
    const c = "ACGT"[Math.floor(next() * 4)];
    if (out.length >= 2 && out[out.length - 1] === out[out.length - 2] && out[out.length - 1] === c) continue;
    out.push(c);
  }
  return out.join("");
}
