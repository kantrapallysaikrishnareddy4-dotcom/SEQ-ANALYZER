"""
matcher.py  -  "Which known protein does my sequence match, and what is different?"

What it does
  1. Takes your sequence (protein, or DNA/RNA which it translates in all 6 reading frames).
  2. Aligns it against the reference proteins stored in reference.db.
  3. Lists the differences from the reference (substitutions, deletions, insertions).
  4. Looks each difference up in the known-variant table (with disease information).

Think of it like this: the database is a library of "correct copies" of genes.
Your sequence is a page someone handed you. We find which book it came from,
then lay the pages side by side and circle every place that differs.

Everything here works with plain text and numbers. It never touches the web page.
"""

import math
import os
import sqlite3
from collections import Counter, defaultdict

import sequence_analyzer as sa      # reuse the codon table from your original script

try:                                 # Biopython is fast. If it is missing we use a slower backup.
    if os.environ.get("SEQ_NO_BIO"):
        raise ImportError("Biopython switched off by SEQ_NO_BIO")
    from Bio.Align import PairwiseAligner, substitution_matrices
    HAVE_BIO = True
except Exception:                    # pragma: no cover
    HAVE_BIO = False

DB_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference.db")

# ---------------- settings ----------------
MIN_SEGMENT = 12          # ignore translated stretches shorter than this (amino acids)
MIN_PROTEIN = 8           # shortest protein we will search
MAX_DNA = 200_000         # longest DNA/RNA we will translate
MAX_SEGMENTS = 40         # align at most this many translated stretches
CANDIDATES = 8            # references to score per stretch after the quick k-mer filter
KMER = 3
SHOW_E = 1e-3             # only show matches with an E-value at least this good
MAX_HITS = 5
CELL_LIMIT = 60_000_000   # skip alignments bigger than this (reference length x query length)
PY_CELL_LIMIT = 4_000_000 # the backup aligner is slower, so its limit is smaller
ALN_ROWS = 40             # alignment rows (60 columns each) sent to the page

# Statistics for BLOSUM62 with gap open 11 / extend 1 (the BLAST protein defaults)
LAMBDA, K = 0.267, 0.041
GAP_OPEN, GAP_EXTEND = 12, 1      # cost of a 1-residue gap is 12, each extra residue costs 1 more

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln", "E": "Glu",
    "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys", "M": "Met", "F": "Phe",
    "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp", "Y": "Tyr", "V": "Val", "X": "Xaa",
}

CODON = {codon.replace("U", "T"): aa for codon, aa in sa.CODON_TABLE.items()}
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


# ---------------------------------------------------------------
# BLOSUM62 (copied from Biopython, so the backup aligner scores the same)
# ---------------------------------------------------------------
_BLOSUM_TEXT = """     A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
A    4  -1  -2  -2   0  -1  -1   0  -2  -1  -1  -1  -1  -2  -1   1   0  -3  -2   0  -2  -1   0  -4
R   -1   5   0  -2  -3   1   0  -2   0  -3  -2   2  -1  -3  -2  -1  -1  -3  -2  -3  -1   0  -1  -4
N   -2   0   6   1  -3   0   0   0   1  -3  -3   0  -2  -3  -2   1   0  -4  -2  -3   3   0  -1  -4
D   -2  -2   1   6  -3   0   2  -1  -1  -3  -4  -1  -3  -3  -1   0  -1  -4  -3  -3   4   1  -1  -4
C    0  -3  -3  -3   9  -3  -4  -3  -3  -1  -1  -3  -1  -2  -3  -1  -1  -2  -2  -1  -3  -3  -2  -4
Q   -1   1   0   0  -3   5   2  -2   0  -3  -2   1   0  -3  -1   0  -1  -2  -1  -2   0   3  -1  -4
E   -1   0   0   2  -4   2   5  -2   0  -3  -3   1  -2  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
G    0  -2   0  -1  -3  -2  -2   6  -2  -4  -4  -2  -3  -3  -2   0  -2  -2  -3  -3  -1  -2  -1  -4
H   -2   0   1  -1  -3   0   0  -2   8  -3  -3  -1  -2  -1  -2  -1  -2  -2   2  -3   0   0  -1  -4
I   -1  -3  -3  -3  -1  -3  -3  -4  -3   4   2  -3   1   0  -3  -2  -1  -3  -1   3  -3  -3  -1  -4
L   -1  -2  -3  -4  -1  -2  -3  -4  -3   2   4  -2   2   0  -3  -2  -1  -2  -1   1  -4  -3  -1  -4
K   -1   2   0  -1  -3   1   1  -2  -1  -3  -2   5  -1  -3  -1   0  -1  -3  -2  -2   0   1  -1  -4
M   -1  -1  -2  -3  -1   0  -2  -3  -2   1   2  -1   5   0  -2  -1  -1  -1  -1   1  -3  -1  -1  -4
F   -2  -3  -3  -3  -2  -3  -3  -3  -1   0   0  -3   0   6  -4  -2  -2   1   3  -1  -3  -3  -1  -4
P   -1  -2  -2  -1  -3  -1  -1  -2  -2  -3  -3  -1  -2  -4   7  -1  -1  -4  -3  -2  -2  -1  -2  -4
S    1  -1   1   0  -1   0   0   0  -1  -2  -2   0  -1  -2  -1   4   1  -3  -2  -2   0   0   0  -4
T    0  -1   0  -1  -1  -1  -1  -2  -2  -1  -1  -1  -1  -2  -1   1   5  -2  -2   0  -1  -1   0  -4
W   -3  -3  -4  -4  -2  -2  -3  -2  -2  -3  -2  -3  -1   1  -4  -3  -2  11   2  -3  -4  -3  -2  -4
Y   -2  -2  -2  -3  -2  -1  -2  -3   2  -1  -1  -2  -1   3  -3  -2  -2   2   7  -1  -3  -2  -1  -4
V    0  -3  -3  -3  -1  -2  -2  -3  -3   3   1  -2   1  -1  -2  -2   0  -3  -1   4  -3  -2  -1  -4
B   -2  -1   3   4  -3   0   1  -1   0  -3  -4   0  -3  -3  -2   0  -1  -4  -3  -3   4   1  -1  -4
Z   -1   0   0   1  -3   3   4  -2   0  -3  -3   1  -1  -3  -1   0  -1  -3  -2  -2   1   4  -1  -4
X    0  -1  -1  -1  -2  -1  -1  -1  -1  -1  -1  -1  -1  -1  -2   0   0  -2  -1  -1  -1  -1  -1  -4
*   -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4  -4   1"""


def _parse_blosum(text):
    rows = [ln.split() for ln in text.strip().splitlines()]
    cols = rows[0]
    table = {}
    for row in rows[1:]:
        letter = row[0]
        for c, v in zip(cols, row[1:]):
            table[(letter, c)] = int(v)
    return table


BLOSUM = _parse_blosum(_BLOSUM_TEXT)


# ---------------------------------------------------------------
# Aligners: a fast one (Biopython) and a slow backup (pure Python)
# ---------------------------------------------------------------
if HAVE_BIO:
    _aligner = PairwiseAligner()
    _aligner.mode = "local"
    _aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    _aligner.open_gap_score = -GAP_OPEN
    _aligner.extend_gap_score = -GAP_EXTEND


def _score_py(a, b):
    """Smith-Waterman score only (affine gaps), keeps just two rows in memory."""
    m = len(b)
    NEG = -10 ** 9
    Hprev = [0] * (m + 1)
    Fprev = [NEG] * (m + 1)
    best = 0
    for ai in a:
        Hcur = [0] * (m + 1)
        Fcur = [NEG] * (m + 1)
        e = NEG
        for j in range(1, m + 1):
            e = max(Hcur[j - 1] - GAP_OPEN, e - GAP_EXTEND)
            f = max(Hprev[j] - GAP_OPEN, Fprev[j] - GAP_EXTEND)
            h = Hprev[j - 1] + BLOSUM[(ai, b[j - 1])]
            if e > h:
                h = e
            if f > h:
                h = f
            if h < 0:
                h = 0
            Hcur[j] = h
            Fcur[j] = f
            if h > best:
                best = h
        Hprev, Fprev = Hcur, Fcur
    return best


def _align_py(a, b):
    """Smith-Waterman with traceback. Returns (score, target_blocks, query_blocks)."""
    n, m = len(a), len(b)
    NEG = -10 ** 9
    H = [[0] * (m + 1) for _ in range(n + 1)]
    E = [[NEG] * (m + 1) for _ in range(n + 1)]
    F = [[NEG] * (m + 1) for _ in range(n + 1)]
    best, bi, bj = 0, 0, 0
    for i in range(1, n + 1):
        ai = a[i - 1]
        Hi, Hp, Ei, Fi, Fp = H[i], H[i - 1], E[i], F[i], F[i - 1]
        for j in range(1, m + 1):
            e = max(Hi[j - 1] - GAP_OPEN, Ei[j - 1] - GAP_EXTEND)
            f = max(Hp[j] - GAP_OPEN, Fp[j] - GAP_EXTEND)
            h = Hp[j - 1] + BLOSUM[(ai, b[j - 1])]
            if e > h:
                h = e
            if f > h:
                h = f
            if h < 0:
                h = 0
            Hi[j], Ei[j], Fi[j] = h, e, f
            if h > best:
                best, bi, bj = h, i, j
    if best == 0:
        return 0, [], []

    i, j, state = bi, bj, "H"
    pairs = []
    while True:
        if state == "H":
            h = H[i][j]
            if h == 0:
                break
            if i > 0 and j > 0 and h == H[i - 1][j - 1] + BLOSUM[(a[i - 1], b[j - 1])]:
                pairs.append((i - 1, j - 1))
                i -= 1
                j -= 1
            elif h == E[i][j]:
                state = "E"
            elif h == F[i][j]:
                state = "F"
            else:
                break
        elif state == "E":              # gap in the reference: uses up a query residue
            state = "H" if E[i][j] == H[i][j - 1] - GAP_OPEN else "E"
            j -= 1
        else:                           # gap in the query: uses up a reference residue
            state = "H" if F[i][j] == H[i - 1][j] - GAP_OPEN else "F"
            i -= 1
    pairs.reverse()

    tb, qb = [], []
    for ti, qi in pairs:
        if tb and tb[-1][1] == ti and qb[-1][1] == qi:
            tb[-1] = (tb[-1][0], ti + 1)
            qb[-1] = (qb[-1][0], qi + 1)
        else:
            tb.append((ti, ti + 1))
            qb.append((qi, qi + 1))
    return best, tb, qb


def quick_score(ref, seg):
    if HAVE_BIO:
        if len(ref) * len(seg) > CELL_LIMIT:
            return 0
        return _aligner.score(ref, seg)
    if len(ref) * len(seg) > PY_CELL_LIMIT:
        return 0
    return _score_py(ref, seg)


def full_align(ref, seg):
    if HAVE_BIO:
        aln = _aligner.align(ref, seg)[0]
        t, q = aln.aligned
        return (float(aln.score),
                [(int(x), int(y)) for x, y in t],
                [(int(x), int(y)) for x, y in q])
    return _align_py(ref, seg)


# ---------------------------------------------------------------
# Database access
# ---------------------------------------------------------------
_cache = {"path": None, "mtime": None}


def db_path():
    return os.environ.get("SEQ_DB", DB_PATH)


def _connect(path):
    conn = sqlite3.connect(path)
    conn.row_factory = sqlite3.Row
    return conn


def db_info(path=None):
    """Is the database built, and what is in it?"""
    path = path or db_path()
    empty = {"ready": False, "proteins": 0, "variants": 0, "meta": {}}
    if not os.path.exists(path):
        return empty
    try:
        conn = _connect(path)
        proteins = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
        variants = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
        meta = {r[0]: r[1] for r in conn.execute("SELECT key, value FROM meta").fetchall()}
        conn.close()
    except sqlite3.Error:
        return empty
    return {"ready": proteins > 0, "proteins": proteins, "variants": variants, "meta": meta}


def _load_refs(path):
    """Read all reference proteins into memory once (and again if the file changes)."""
    mtime = os.path.getmtime(path)
    if _cache["path"] == path and _cache["mtime"] == mtime:
        return _cache
    conn = _connect(path)
    refs, index = [], defaultdict(set)
    for i, row in enumerate(conn.execute(
            "SELECT accession, gene, name, organism, length, sequence FROM proteins")):
        seq = row["sequence"]
        refs.append({"accession": row["accession"], "gene": row["gene"], "name": row["name"],
                     "organism": row["organism"], "length": len(seq), "sequence": seq})
        for k in range(len(seq) - KMER + 1):
            index[seq[k:k + KMER]].add(i)
    conn.close()
    _cache.update({"path": path, "mtime": mtime, "refs": refs, "index": index,
                   "total": sum(r["length"] for r in refs)})
    return _cache


# ---------------------------------------------------------------
# Turning DNA into protein pieces
# ---------------------------------------------------------------
def translate_segments(dna):
    """Translate all 6 frames and cut each at stop codons. Returns the pieces worth searching."""
    L = len(dna)
    out = []
    for strand, s in (("+", dna), ("-", revcomp(dna))):
        for frame in range(3):
            n = (L - frame) // 3
            prot = "".join(CODON[s[frame + 3 * i: frame + 3 * i + 3]] for i in range(n))
            pos = 0
            for piece in prot.split("*"):
                if len(piece) >= MIN_SEGMENT:
                    out.append({"seq": piece, "strand": strand, "frame": frame,
                                "aa_start": pos, "stop_after": pos + len(piece) < n, "L": L})
                pos += len(piece) + 1
    return out


def nt_range(L, strand, frame, aa_i, aa_j):
    """Amino-acid indexes [aa_i, aa_j) of a frame -> nucleotide range on the ORIGINAL sequence."""
    if strand == "+":
        return frame + aa_i * 3 + 1, frame + aa_j * 3
    return L - (frame + aa_j * 3) + 1, L - (frame + aa_i * 3)


def candidates(seg, data):
    """Quick filter: which references share the most 3-letter words with this piece?"""
    counts = Counter()
    for k in {seg[i:i + KMER] for i in range(len(seg) - KMER + 1)}:
        for ri in data["index"].get(k, ()):
            counts[ri] += 1
    return [ri for ri, _ in counts.most_common(CANDIDATES)]


# ---------------------------------------------------------------
# Reading the alignment: differences, names, text view
# ---------------------------------------------------------------
def call_differences(ref, seg, tb, qb):
    diffs = []
    for idx, ((ts, te), (qs, qe)) in enumerate(zip(tb, qb)):
        if idx > 0:
            pte, pqe = tb[idx - 1][1], qb[idx - 1][1]
            gt, gq = ts - pte, qs - pqe
            if gt > 0 and gq > 0:
                diffs.append({"type": "delins", "start": pte + 1, "end": ts,
                              "ref": ref[pte:ts], "alt": seg[pqe:qs]})
            elif gt > 0:
                diffs.append({"type": "deletion", "start": pte + 1, "end": ts,
                              "ref": ref[pte:ts], "alt": ""})
            elif gq > 0:
                diffs.append({"type": "insertion", "start": pte, "end": pte + 1,
                              "ref": "", "alt": seg[pqe:qs]})
        for k in range(te - ts):
            r, q = ref[ts + k], seg[qs + k]
            if r != q:
                diffs.append({"type": "substitution", "start": ts + k + 1, "end": ts + k + 1,
                              "ref": r, "alt": q})
    for d in diffs:
        normalize_indel(ref, d)
    diffs.sort(key=lambda d: (d["start"], d["end"]))
    return diffs


def normalize_indel(ref, d):
    """If a deletion/insertion sits inside a repeat (like LL or VV), the alignment could place it
    in several equivalent spots. The standard convention is to report the most C-terminal one."""
    if d["type"] == "deletion":
        while d["end"] < len(ref) and ref[d["start"] - 1] == ref[d["end"]]:
            d["start"] += 1
            d["end"] += 1
        d["ref"] = ref[d["start"] - 1:d["end"]]
    elif d["type"] == "insertion":
        p, alt = d["start"], d["alt"]
        while p < len(ref) and ref[p] == alt[0]:
            alt = alt[1:] + alt[0]
            p += 1
        d["start"], d["end"], d["alt"] = p, p + 1, alt


def _three(letters):
    return "".join(AA3.get(c, "Xaa") for c in letters)


def hgvs_name(d, ref):
    """Standard protein change names, e.g. p.Glu7Val, p.Phe508del."""
    t = d["type"]
    if t == "substitution":
        return f"p.{_three(d['ref'])}{d['start']}{_three(d['alt'])}"
    if t == "deletion":
        first = f"{_three(ref[d['start'] - 1])}{d['start']}"
        if d["start"] == d["end"]:
            return f"p.{first}del"
        return f"p.{first}_{_three(ref[d['end'] - 1])}{d['end']}del"
    if t == "insertion":
        p = d["start"]
        left = f"{_three(ref[p - 1])}{p}"
        if p < len(ref):
            return f"p.{left}_{_three(ref[p])}{p + 1}ins{_three(d['alt'])}"
        return f"p.{left}ins{_three(d['alt'])}"
    first = f"{_three(ref[d['start'] - 1])}{d['start']}"
    if d["start"] == d["end"]:
        return f"p.{first}delins{_three(d['alt'])}"
    return f"p.{first}_{_three(ref[d['end'] - 1])}{d['end']}delins{_three(d['alt'])}"


EFFECT = {"substitution": "missense (one amino acid changed)",
          "deletion": "deletion", "insertion": "insertion", "delins": "deletion + insertion"}


def alignment_rows(ref, seg, tb, qb, width=60):
    """Build the side-by-side text view: reference on top, query below, | for identical."""
    rcols, qcols = [], []
    for idx, ((ts, te), (qs, qe)) in enumerate(zip(tb, qb)):
        if idx > 0:
            pte, pqe = tb[idx - 1][1], qb[idx - 1][1]
            gt, gq = ts - pte, qs - pqe
            rcols += list(ref[pte:ts]) + ["-"] * gq
            qcols += ["-"] * gt + list(seg[pqe:qs])
        rcols += list(ref[ts:te])
        qcols += list(seg[qs:qe])

    mid = []
    for r, q in zip(rcols, qcols):
        if r == "-" or q == "-":
            mid.append(" ")
        elif r == q:
            mid.append("|")
        elif BLOSUM.get((r, q), -1) > 0:
            mid.append("+")
        else:
            mid.append(" ")

    rows = []
    rpos, qpos = tb[0][0] + 1, qb[0][0] + 1
    for start in range(0, len(rcols), width):
        rpart = "".join(rcols[start:start + width])
        qpart = "".join(qcols[start:start + width])
        rows.append({"ref": rpart, "mid": "".join(mid[start:start + width]), "query": qpart,
                     "ref_start": rpos, "query_start": qpos})
        rpos += sum(1 for c in rpart if c != "-")
        qpos += sum(1 for c in qpart if c != "-")
    return rows


def _split_pipe(text):
    return [p for p in (text or "").split("|") if p]


def _variant_dict(row):
    return {"description": row["description"], "diseases": _split_pipe(row["diseases"]),
            "dbsnp": row["dbsnp"], "feature_id": row["feature_id"]}


def annotate_differences(conn, acc, ref, diffs):
    """Add the HGVS name and any known-variant information to each difference."""
    if not diffs:
        return diffs
    lo = min(d["start"] for d in diffs)
    hi = max(d["end"] for d in diffs)
    by_pos = defaultdict(list)
    for row in conn.execute(
            "SELECT position, ref, alt, description, diseases, dbsnp, feature_id "
            "FROM variants WHERE accession=? AND position BETWEEN ? AND ?", (acc, lo, hi)):
        by_pos[row["position"]].append(row)

    for d in diffs:
        d["hgvs"] = hgvs_name(d, ref)
        d["effect"] = EFFECT[d["type"]]
        rows = by_pos.get(d["start"], [])
        d["known"] = [_variant_dict(r) for r in rows if r["ref"] == d["ref"] and r["alt"] == d["alt"]]
        d["same_site"] = []
        if d["type"] == "substitution":
            d["same_site"] = [
                {"change": f"{_three(r['ref'])}{r['position']}{_three(r['alt'])}", **_variant_dict(r)}
                for r in rows
                if len(r["ref"]) == 1 and len(r["alt"]) == 1 and r["alt"] != d["alt"]][:4]
    return diffs


# ---------------------------------------------------------------
# The main function
# ---------------------------------------------------------------
def match(seq, seq_type, path=None):
    path = path or db_path()
    info = db_info(path)
    if not info["ready"]:
        return {"ready": False}

    data = _load_refs(path)
    notes = []

    if seq_type == "PROTEIN":
        if len(seq) < MIN_PROTEIN:
            return {"ready": True, "hits": [], "refs": len(data["refs"]), "meta": info["meta"],
                    "notes": [f"This sequence is too short to search reliably "
                              f"(needs at least {MIN_PROTEIN} amino acids)."]}
        segments = [{"seq": seq, "strand": "", "frame": 0, "aa_start": 0,
                     "stop_after": False, "L": 0}]
    else:
        dna = seq.replace("U", "T")
        if len(dna) > MAX_DNA:
            return {"ready": True, "error": f"This sequence is too long for the built-in matcher "
                                            f"(limit {MAX_DNA:,} letters)."}
        segments = translate_segments(dna)
        if not segments:
            return {"ready": True, "hits": [], "refs": len(data["refs"]), "meta": info["meta"],
                    "notes": [f"No stretch of at least {MIN_SEGMENT} amino acids without a stop codon "
                              f"was found in any of the 6 reading frames."]}
        if len(segments) > MAX_SEGMENTS:
            segments.sort(key=lambda s: -len(s["seq"]))
            segments = segments[:MAX_SEGMENTS]
            notes.append(f"Many translated stretches were found. Only the {MAX_SEGMENTS} longest were searched.")

    # Step 1: quick scores. Keep the best-scoring piece for each reference protein.
    best = {}                                  # ref index -> (score, segment index)
    for si, seg in enumerate(segments):
        for ri in candidates(seg["seq"], data):
            score = quick_score(data["refs"][ri]["sequence"], seg["seq"])
            if score > 0 and (ri not in best or score > best[ri][0]):
                best[ri] = (score, si)

    # Step 2: significance. E-value = how many matches this good we'd expect by pure chance.
    total = data["total"]
    scored = []
    for ri, (score, si) in best.items():
        m = len(segments[si]["seq"])
        evalue = K * m * total * math.exp(-LAMBDA * score)
        if evalue <= SHOW_E:
            scored.append((score, evalue, ri, si))
    scored.sort(key=lambda t: (-t[0], t[1]))

    # Step 3: full alignment, differences and database lookups for the top hits.
    conn = _connect(path)
    hits = []
    for score, evalue, ri, si in scored[:MAX_HITS]:
        ref = data["refs"][ri]
        seg = segments[si]
        _, tb, qb = full_align(ref["sequence"], seg["seq"])
        if not tb:
            continue

        matches = sum(1 for (ts, te), (qs, _) in zip(tb, qb)
                      for k in range(te - ts) if ref["sequence"][ts + k] == seg["seq"][qs + k])
        aligned_ref = sum(te - ts for ts, te in tb)
        aligned_query = sum(qe - qs for qs, qe in qb)
        gaps = sum(max(tb[i][0] - tb[i - 1][1], 0) + max(qb[i][0] - qb[i - 1][1], 0)
                   for i in range(1, len(tb)))
        aln_len = aligned_ref + sum(max(qb[i][0] - qb[i - 1][1], 0) for i in range(1, len(tb))) \
            + sum(max(tb[i][0] - tb[i - 1][1], 0) for i in range(1, len(tb)))
        ref_start, ref_end = tb[0][0] + 1, tb[-1][1]
        q_start, q_end = qb[0][0] + 1, qb[-1][1]

        diffs = call_differences(ref["sequence"], seg["seq"], tb, qb)
        annotate_differences(conn, ref["accession"], ref["sequence"], diffs)

        flags = []
        if seg["stop_after"] and len(seg["seq"]) - q_end <= 3 and ref["length"] - ref_end >= 5:
            nxt = ref["sequence"][ref_end]
            flags.append({
                "kind": "truncated",
                "text": (f"Your sequence reaches a stop codon right after reference residue {ref_end}, "
                         f"but the reference protein continues to residue {ref['length']}. "
                         f"If this is real, the protein would be cut short "
                         f"(a nonsense change, likely p.{_three(nxt)}{ref_end + 1}*). "
                         f"It could also be a frameshift or an incomplete sequence."),
            })
        if ref_start > 1:
            flags.append({"kind": "partial", "text": f"Your sequence starts at reference residue {ref_start}."})
        if ref_end < ref["length"] and not flags_has(flags, "truncated"):
            flags.append({"kind": "partial", "text": f"Your sequence ends at reference residue {ref_end} of {ref['length']}."})

        rows = alignment_rows(ref["sequence"], seg["seq"], tb, qb)
        hit = {
            "accession": ref["accession"], "gene": ref["gene"], "name": ref["name"],
            "organism": ref["organism"], "ref_length": ref["length"],
            "score": round(score, 1), "bits": round((LAMBDA * score - math.log(K)) / math.log(2), 1),
            "evalue": evalue,
            "identity": round(matches / aln_len * 100, 1) if aln_len else 0,
            "matches": matches, "aligned_len": aln_len, "gaps": gaps,
            "query_cov": round(aligned_query / len(seg["seq"]) * 100, 1),
            "query_len": len(seg["seq"]),
            "ref_start": ref_start, "ref_end": ref_end, "query_start": q_start, "query_end": q_end,
            "strand": seg["strand"], "frame": seg["frame"] + 1 if seg["strand"] else None,
            "nt_start": None, "nt_end": None,
            "flags": flags, "differences": diffs,
            "alignment": rows[:ALN_ROWS], "alignment_truncated": len(rows) > ALN_ROWS,
        }
        if seg["strand"]:
            hit["nt_start"], hit["nt_end"] = nt_range(
                seg["L"], seg["strand"], seg["frame"],
                seg["aa_start"] + qb[0][0], seg["aa_start"] + qb[-1][1])

        hit["diseases"] = [
            {"acronym": r["acronym"], "name": r["name"], "mim": r["mim"]}
            for r in conn.execute(
                "SELECT acronym, name, mim FROM protein_diseases WHERE accession=?", (ref["accession"],))]
        counts = conn.execute(
            "SELECT COUNT(*), SUM(CASE WHEN diseases != '' THEN 1 ELSE 0 END) "
            "FROM variants WHERE accession=?", (ref["accession"],)).fetchone()
        hit["variants_on_file"] = counts[0] or 0
        hit["disease_variants_on_file"] = counts[1] or 0
        hits.append(hit)
    conn.close()

    if not hits:
        notes.append("No significant match in this database. It only holds the proteins listed in "
                     "the database info, so a sequence that is not among them cannot match.")
    return {"ready": True, "hits": hits, "refs": len(data["refs"]), "meta": info["meta"],
            "notes": notes, "engine": "biopython" if HAVE_BIO else "python-backup"}


def flags_has(flags, kind):
    return any(f["kind"] == kind for f in flags)
