"""
app.py  -  connects your Python analyzer to the webpage

How the pieces talk to each other:

  Browser (HTML / CSS / JS)  --- sends the sequence --->  Flask (this file)
                                                             |
                                                             |  uses your original code
                                                             v
                                                     sequence_analyzer.py
  Browser  <--- gets the results back (as JSON) ---  Flask

Run it:   python app.py      then open   http://127.0.0.1:5000
"""

from collections import Counter

from flask import Flask, jsonify, request, send_from_directory

import sequence_analyzer as sa      # <-- YOUR original script, reused here
import matcher                      # <-- finds database matches and mutations

# static_url_path="" lets the page load style.css and script.js directly
app = Flask(__name__, static_folder="static", static_url_path="")

MAX_INPUT = 2_000_000               # refuse absurdly large pastes

# Your original codon table uses U (RNA). The page works with T (DNA), so convert it once.
CODON_TABLE = {codon.replace("U", "T"): aa for codon, aa in sa.CODON_TABLE.items()}

# Kyte-Doolittle scale (water-hating > 0, water-loving < 0), used for GRAVY
HYDROPATHY = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
PROTEIN_LETTERS = "ACDEFGHIKLMNPQRSTVWY"
FRAME_CHOICES = {"+1", "+2", "+3", "-1", "-2", "-3"}


# ---------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------

def parse_input(raw):
    """Pull out the FASTA name and the first sequence. Returns (name, letters, had_extra_records)."""
    name = ""
    kept = []
    started = False
    extra = False
    for line in raw.splitlines():
        if line.strip().startswith(">"):
            if started:             # a second FASTA entry: ignore it
                extra = True
                break
            name = line.strip()[1:].strip()
            started = True
        else:
            kept.append(line)
            if line.strip():
                started = True
    # your original clean_sequence keeps only the letters
    return name, sa.clean_sequence("\n".join(kept)), extra


def translate(codons):
    return "".join(CODON_TABLE[c] for c in codons)


def find_orfs(dna):
    """Every ORF (ATG ... stop) in all 3 frames on both strands."""
    length = len(dna)
    found = []
    strands = [("+", dna), ("-", sa.reverse_complement(dna, "DNA"))]

    for sign, strand in strands:
        for frame in range(3):
            codons = sa.split_into_codons(strand[frame:])
            i = 0
            while i < len(codons):
                if codons[i] != "ATG":
                    i += 1
                    continue
                j = i
                protein = ""
                while j < len(codons) and CODON_TABLE[codons[j]] != "*":
                    protein += CODON_TABLE[codons[j]]
                    j += 1
                if j >= len(codons):        # no stop codon before the end: not a full ORF
                    break
                start = frame + i * 3 + 1
                end = frame + (j + 1) * 3
                if sign == "-":             # convert to positions on the original sequence
                    start, end = length - end + 1, length - start + 1
                found.append({"strand": sign, "frame": frame + 1, "start": start,
                              "end": end, "aa": len(protein), "protein": protein})
                i = j + 1

    found.sort(key=lambda o: (-o["aa"], o["start"]))
    return found


# ---------------------------------------------------------------
# The two analyses (they return plain dictionaries -> JSON)
# ---------------------------------------------------------------

def analyze_nucleic(seq, seq_type, frame_choice, min_orf):
    is_rna = seq_type == "RNA"
    dna = seq.replace("U", "T") if is_rna else seq    # work in T-letters internally
    length = len(dna)

    counts = Counter(dna)
    a, c, g, t = (counts.get(b, 0) for b in "ACGT")

    def pct(n):
        return n / length * 100

    sign = frame_choice[0]
    frame = int(frame_choice[1]) - 1
    strand = dna if sign == "+" else sa.reverse_complement(dna, "DNA")
    codons = sa.split_into_codons(strand[frame:])

    starts, stops = [], []
    for i, codon in enumerate(codons):
        info = {"codon": i + 1, "nt": frame + i * 3 + 1}
        if codon == "ATG":
            starts.append(info)
        if CODON_TABLE[codon] == "*":
            stops.append({**info, "text": codon})

    usage = sorted(Counter(codons).items(), key=lambda kv: (-kv[1], kv[0]))

    return {
        "kind": "nucleic", "type": seq_type, "isRNA": is_rna, "L": length,
        "counts": {"A": a, "C": c, "G": g, "T": t},
        "pct": {"A": pct(a), "C": pct(c), "G": pct(g), "T": pct(t)},
        "gc": sa.gc_content(dna),                      # <-- from your original script
        "at": pct(a + t),
        "purines": a + g, "pyrimidines": c + t,
        "dna": dna,
        "comp": sa.complement(dna, "DNA"),             # <-- from your original script
        "rc": sa.reverse_complement(dna, "DNA"),       # <-- from your original script
        "frameChoice": frame_choice, "strandSign": sign, "frame": frame,
        "codons": codons, "leftover": len(strand) - frame - len(codons) * 3,
        "usage": usage, "starts": starts, "stops": stops,
        "protein": translate(codons),
        "orfs": [o for o in find_orfs(dna) if o["aa"] >= min_orf],
    }


def analyze_protein(seq):
    length = len(seq)
    counts = Counter(seq)
    groups = [{"name": name, "n": n, "pct": n / length * 100}
              for name, n in sa.protein_groups(seq).items()]   # <-- from your original script
    return {
        "kind": "protein", "L": length, "counts": dict(counts),
        "mass": sa.protein_mass(seq),                          # <-- from your original script
        "groups": groups,
        "net": counts["K"] + counts["R"] - counts["D"] - counts["E"],
        "gravy": sum(HYDROPATHY[aa] for aa in seq) / length,
        "table": sorted(counts.items(), key=lambda kv: (-kv[1], kv[0])),
    }


# ---------------------------------------------------------------
# Web routes: the "addresses" the browser can talk to
# ---------------------------------------------------------------

@app.route("/")
def home():
    return send_from_directory(app.static_folder, "index.html")


@app.post("/analyze")
def analyze():
    data = request.get_json(silent=True) or {}
    raw = str(data.get("sequence", ""))
    if len(raw) > MAX_INPUT:
        return jsonify({"error": "That input is too large (limit is 2 million characters)."})

    name, seq, extra = parse_input(raw)
    if not seq:
        return jsonify({"empty": True})

    chosen = data.get("type", "auto")
    chosen_by_user = chosen in ("DNA", "RNA", "PROTEIN")
    seq_type = chosen if chosen_by_user else sa.detect_type(seq)

    frame_choice = data.get("frame", "+1")
    if frame_choice not in FRAME_CHOICES:
        frame_choice = "+1"
    try:
        min_orf = max(1, int(data.get("minOrf", 1)))
    except (TypeError, ValueError):
        min_orf = 1

    # Make sure the letters really fit the type
    if seq_type == "UNKNOWN":
        bad = ", ".join(sorted(set(seq) - set(PROTEIN_LETTERS + "U")))
        return jsonify({"error": f"Could not recognise this sequence. Unexpected letters: {bad}. "
                                 "If it mixes T and U, it is not a single DNA or RNA sequence."})
    allowed = {"DNA": "ACGT", "RNA": "ACGU", "PROTEIN": PROTEIN_LETTERS}[seq_type]
    bad = sorted(set(seq) - set(allowed))
    if bad:
        return jsonify({"error": f"This does not look like {seq_type}. These letters don't belong: "
                                 f"{', '.join(bad)}. Check the Type menu or your sequence."})

    notes = []
    if extra:
        notes.append({"level": "warn",
                      "text": "Your input has more than one FASTA entry. Only the first one was analyzed."})
    if not chosen_by_user and seq_type != "PROTEIN" and set(seq) <= set("ACGT") and len(seq) < 30:
        notes.append({"level": "",
                      "text": f"Detected as {seq_type}. If this is really a protein made only of "
                              "A, C, G, T, pick Protein in the Type menu."})

    if seq_type == "PROTEIN":
        result = analyze_protein(seq)
    else:
        result = analyze_nucleic(seq, seq_type, frame_choice, min_orf)

    return jsonify({"name": name, "notes": notes, "minOrf": min_orf, "result": result})


@app.get("/db-info")
def db_info():
    """Tells the page whether the reference database has been built."""
    return jsonify(matcher.db_info())


@app.post("/match")
def match():
    """Align the sequence against the reference database and list the differences."""
    data = request.get_json(silent=True) or {}
    raw = str(data.get("sequence", ""))
    if len(raw) > MAX_INPUT:
        return jsonify({"error": "That input is too large."})

    _, seq, _ = parse_input(raw)
    if not seq:
        return jsonify({"empty": True})

    chosen = data.get("type", "auto")
    seq_type = chosen if chosen in ("DNA", "RNA", "PROTEIN") else sa.detect_type(seq)
    allowed = {"DNA": "ACGT", "RNA": "ACGU", "PROTEIN": PROTEIN_LETTERS}.get(seq_type)
    if allowed is None or set(seq) - set(allowed):
        return jsonify({"error": "Could not search: the sequence letters do not fit the type."})

    return jsonify(matcher.match(seq, seq_type))


if __name__ == "__main__":
    # debug=True reloads the server automatically when you save changes
    app.run(debug=True)
