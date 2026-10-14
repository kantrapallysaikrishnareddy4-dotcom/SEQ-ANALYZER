"""
Sequence Analyzer: DNA / RNA / Protein
--------------------------------------
Paste a sequence and the program will:
  - detect if it is DNA, RNA or protein
  - count every nucleotide (A, T/U, G, C) or amino acid
  - calculate GC content and AT content
  - split DNA/RNA into codons and translate them to protein
  - find start (ATG/AUG) and stop codons
  - give complement, reverse complement and transcription (for DNA)
"""

from collections import Counter

# ---------------------------------------------------------------
# 1. REFERENCE DATA (like a dictionary book the program looks up)
# ---------------------------------------------------------------

# Codon table: each 3-letter RNA codon -> amino acid (1-letter code)
# "*" means STOP
CODON_TABLE = {
    "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
    "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
    "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
    "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
    "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
    "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "UAU": "Y", "UAC": "Y", "UAA": "*", "UAG": "*",
    "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
    "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

AMINO_ACID_NAMES = {
    "A": "Alanine", "R": "Arginine", "N": "Asparagine", "D": "Aspartic acid",
    "C": "Cysteine", "E": "Glutamic acid", "Q": "Glutamine", "G": "Glycine",
    "H": "Histidine", "I": "Isoleucine", "L": "Leucine", "K": "Lysine",
    "M": "Methionine", "F": "Phenylalanine", "P": "Proline", "S": "Serine",
    "T": "Threonine", "W": "Tryptophan", "Y": "Tyrosine", "V": "Valine",
    "*": "STOP",
}

# Average mass of each amino acid inside a protein chain (Daltons)
RESIDUE_MASS = {
    "A": 71.0788, "R": 156.1875, "N": 114.1038, "D": 115.0886,
    "C": 103.1388, "E": 129.1155, "Q": 128.1307, "G": 57.0519,
    "H": 137.1411, "I": 113.1594, "L": 113.1594, "K": 128.1741,
    "M": 131.1926, "F": 147.1766, "P": 97.1167, "S": 87.0782,
    "T": 101.1051, "W": 186.2132, "Y": 163.1760, "V": 99.1326,
}
WATER = 18.01524

DNA_COMPLEMENT = {"A": "T", "T": "A", "G": "C", "C": "G"}
RNA_COMPLEMENT = {"A": "U", "U": "A", "G": "C", "C": "G"}


# ---------------------------------------------------------------
# 2. CLEANING + DETECTING
# ---------------------------------------------------------------

def clean_sequence(raw_text):
    """Remove FASTA header lines (starting with >), spaces, numbers, newlines."""
    lines = [ln for ln in raw_text.splitlines() if not ln.startswith(">")]
    text = "".join(lines).upper()
    return "".join(ch for ch in text if ch.isalpha())


def detect_type(seq):
    """Guess whether the sequence is DNA, RNA or protein."""
    letters = set(seq)
    if letters <= set("ACGT"):
        return "DNA"
    if letters <= set("ACGU"):
        return "RNA"
    if letters <= set("ACDEFGHIKLMNPQRSTVWY"):
        return "PROTEIN"
    return "UNKNOWN"


# ---------------------------------------------------------------
# 3. DNA / RNA TOOLS
# ---------------------------------------------------------------

def nucleotide_counts(seq, seq_type):
    """Count each nucleotide and return counts + percentages."""
    bases = "ACGT" if seq_type == "DNA" else "ACGU"
    counts = Counter(seq)
    total = len(seq)
    return {b: (counts[b], counts[b] / total * 100) for b in bases}


def gc_content(seq):
    """GC content = (G + C) / total length x 100"""
    return (seq.count("G") + seq.count("C")) / len(seq) * 100


def to_rna(dna):
    """DNA -> RNA: replace T with U."""
    return dna.replace("T", "U")


def complement(seq, seq_type):
    table = DNA_COMPLEMENT if seq_type == "DNA" else RNA_COMPLEMENT
    return "".join(table[b] for b in seq)


def reverse_complement(seq, seq_type):
    return complement(seq, seq_type)[::-1]


def split_into_codons(rna):
    """Cut the sequence into groups of 3 letters (leftover letters are dropped)."""
    return [rna[i:i + 3] for i in range(0, len(rna) - 2, 3)]


def translate(codons):
    """Turn each codon into its amino acid."""
    return "".join(CODON_TABLE[c] for c in codons)


def find_positions(codons, targets):
    """Return codon positions (1-based) where a codon matches one of the targets."""
    return [(i + 1, c) for i, c in enumerate(codons) if c in targets]


def find_first_orf(rna):
    """Find the first ORF: starts at AUG, ends at a stop codon. Checks all 3 frames."""
    for frame in range(3):
        codons = split_into_codons(rna[frame:])
        for i, c in enumerate(codons):
            if c == "AUG":
                protein = ""
                for j in range(i, len(codons)):
                    aa = CODON_TABLE[codons[j]]
                    if aa == "*":
                        return frame + 1, i * 3 + frame + 1, protein
                    protein += aa
                # reached the end with no stop codon, keep looking
    return None


# ---------------------------------------------------------------
# 4. PROTEIN TOOLS
# ---------------------------------------------------------------

def protein_mass(protein):
    return sum(RESIDUE_MASS[a] for a in protein) + WATER


def protein_groups(protein):
    groups = {
        "Hydrophobic (A,V,I,L,M,F,W,P)": "AVILMFWP",
        "Polar (S,T,N,Q,C,Y,G)": "STNQCYG",
        "Positive charge (K,R,H)": "KRH",
        "Negative charge (D,E)": "DE",
    }
    return {name: sum(protein.count(a) for a in letters)
            for name, letters in groups.items()}


# ---------------------------------------------------------------
# 5. REPORTS (what gets printed)
# ---------------------------------------------------------------

def line(title=""):
    print("\n" + "=" * 55)
    if title:
        print(title)
        print("=" * 55)


def report_nucleic(seq, seq_type):
    line(f"{seq_type} ANALYSIS")
    print(f"Length: {len(seq)} nucleotides")

    print("\n-- Nucleotide counts --")
    for base, (n, pct) in nucleotide_counts(seq, seq_type).items():
        print(f"  {base}: {n:>6}   ({pct:.2f}%)")

    gc = gc_content(seq)
    print(f"\nGC content: {gc:.2f}%")
    print(f"AT/AU content: {100 - gc:.2f}%")

    # Work with RNA from here on (codon table uses U)
    rna = to_rna(seq) if seq_type == "DNA" else seq

    if seq_type == "DNA":
        print("\n-- DNA conversions --")
        print(f"  Complement:         {complement(seq, 'DNA')}")
        print(f"  Reverse complement: {reverse_complement(seq, 'DNA')}")
        print(f"  Transcribed RNA:    {rna}")
    else:
        print("\n-- RNA conversions --")
        print(f"  Complement:         {complement(seq, 'RNA')}")
        print(f"  Reverse complement: {reverse_complement(seq, 'RNA')}")

    codons = split_into_codons(rna)
    print(f"\n-- Codons (reading frame 1) --")
    print(f"  Total codons: {len(codons)}")
    leftover = len(rna) % 3
    if leftover:
        print(f"  (last {leftover} letter(s) ignored, not a full codon)")
    print("  Codon list: " + " ".join(codons))

    print("\n-- Codon usage (how often each codon appears) --")
    for codon, n in sorted(Counter(codons).items()):
        print(f"  {codon} -> {CODON_TABLE[codon]} ({AMINO_ACID_NAMES[CODON_TABLE[codon]]}): {n}")

    print("\n-- Start / stop codons --")
    starts = find_positions(codons, {"AUG"})
    stops = find_positions(codons, {"UAA", "UAG", "UGA"})
    print(f"  Start codons (AUG): {len(starts)} at codon position(s) {[p for p, _ in starts]}")
    print(f"  Stop codons:        {len(stops)} -> {stops}")

    print("\n-- Translation to protein (frame 1) --")
    protein = translate(codons)
    print(f"  {protein}   (* = stop)")

    orf = find_first_orf(rna)
    print("\n-- First ORF (starts at AUG, ends at stop) --")
    if orf:
        frame, start_pos, orf_protein = orf
        print(f"  Frame {frame}, starts at nucleotide {start_pos}")
        print(f"  Protein: {orf_protein}  ({len(orf_protein)} amino acids)")
    else:
        print("  No complete ORF found.")


def report_protein(seq):
    line("PROTEIN ANALYSIS")
    print(f"Length: {len(seq)} amino acids")

    print("\n-- Amino acid counts --")
    counts = Counter(seq)
    for aa, n in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {aa} {AMINO_ACID_NAMES[aa]:<14} {n:>4}  ({n / len(seq) * 100:.2f}%)")

    print(f"\nApproximate molecular weight: {protein_mass(seq):,.2f} Da "
          f"({protein_mass(seq) / 1000:.2f} kDa)")

    print("\n-- Amino acid groups --")
    for name, n in protein_groups(seq).items():
        print(f"  {name}: {n}  ({n / len(seq) * 100:.1f}%)")

    net = seq.count("K") + seq.count("R") - seq.count("D") - seq.count("E")
    print(f"\nRough net charge (K+R minus D+E): {net:+d}")


# ---------------------------------------------------------------
# 6. MAIN PROGRAM
# ---------------------------------------------------------------

def analyze(raw_text):
    seq = clean_sequence(raw_text)
    if not seq:
        print("No sequence found. Please paste letters like ATGGCC...")
        return

    seq_type = detect_type(seq)
    print(f"\nDetected type: {seq_type}")

    if seq_type in ("DNA", "RNA"):
        report_nucleic(seq, seq_type)
    elif seq_type == "PROTEIN":
        report_protein(seq)
    else:
        bad = sorted(set(seq) - set("ACDEFGHIKLMNPQRSTVWYU"))
        print(f"Could not recognise this sequence. Unexpected letters: {bad}")


def main():
    print("DNA / RNA / Protein Analyzer")
    print("Paste your sequence (FASTA is fine). Press Enter twice when done.")
    print("Type 'demo' to try an example, or 'quit' to exit.\n")

    while True:
        lines = []
        while True:
            text = input("> " if not lines else "")
            if text.strip().lower() == "quit":
                return
            if text.strip().lower() == "demo" and not lines:
                lines = ["ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"]
                break
            if text == "":
                break
            lines.append(text)

        if lines:
            analyze("\n".join(lines))
        print("\nPaste another sequence, or type 'quit' to exit.\n")


if __name__ == "__main__":
    main()
