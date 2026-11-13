"""Model of mfeprimer's hairpin structure search (black-box reconstruction, work in progress)."""
C = {"A": "T", "T": "A", "G": "C", "C": "G"}


def comp(a, b):
    return C.get(a) == b


def walk(s, K, minl=3, maxgap=2):
    """Walk diagonal K (column c pairs s[c] with s[K-c]) from the outer end inward; returns the final segment."""
    n = len(s)
    c = max(0, K - n + 1)
    start, cols = c, []
    while (K - c) - c - 1 >= minl:
        ok = comp(s[c], s[K - c])
        cols.append(ok)
        if cols.count(False) > maxgap:
            start, cols = c + 1, []
        c += 1
    while cols and not cols[-1]:
        cols.pop()
    if not cols or not any(cols):
        return None
    pairs = [(start + k, K - start - k) for k, ok in enumerate(cols) if ok]
    return {"K": K, "start": start, "cols": cols, "pairs": pairs, "score": cols.count(True) - cols.count(False)}


def struct(n, pairs):
    st = ["-"] * n
    for i, j in pairs:
        st[i], st[j] = "/", "\\"
    return "".join(st)


def walk2(s, K, minl=3, maxgap=2, mindh=2):
    """H3': inward walk; on gap overflow the segment so far is emitted trimmed at both ends; the segment that
    reaches the centre is emitted trimmed only at the loop side. Returns all candidates."""
    n = len(s)
    c = max(0, K - n + 1)
    start, cols, out = c, [], []

    def emit(st, cs, lead, at):
        cs = list(cs)
        while cs and not cs[-1]:
            cs.pop()
        if not lead:   # a loop-side piece (gaps + helix) that does not raise the score joins the loop
            while cs and False in cs:
                k = len(cs)
                while k and cs[k - 1]:
                    k -= 1
                g = k
                while g and not cs[g - 1]:
                    g -= 1
                if (len(cs) - k) - (k - g) > 0:
                    break
                del cs[g:]
        if lead:
            while cs and not cs[0]:
                cs.pop(0)
                st += 1
        if any(cs):
            pairs = [(st + k, K - st - k) for k, ok in enumerate(cs) if ok]
            sc = cs.count(True) - (0 if lead else cs.count(False))
            out.append({"K": K, "pairs": pairs, "score": sc, "final": not lead, "at": at})

    while (K - c) - c - 1 >= minl:
        cols.append(comp(s[c], s[K - c]))
        if cols.count(False) > maxgap:
            emit(start, cols[:-1], True, c)
            start, cols = c + 1, []
        c += 1
    emit(start, cols, False, c - 1)
    return out


def search(s, minl=3, maxgap=2, mindh=2, pick="dg"):
    """All diagonals; returns the chosen candidate with energy."""
    import hp_energy as E
    best = None
    # diagonals in the order mfeprimer meets them: by centre K//2, the odd one (even loop) first
    for K in sorted(range(2 * len(s) - 1), key=lambda k: (k // 2, -(k % 2))):
        for c in walk2(s, K, minl, maxgap, mindh):
            g, tm = E.energy(s, c["pairs"])
            c.update(dg=g, tm=tm)
            key = {"dg": (g,), "score": (-c["score"],), "scoredg": (-c["score"], g)}[pick]
            if best is None or key < best[0]:
                best = (key, c)
    return best and best[1]
