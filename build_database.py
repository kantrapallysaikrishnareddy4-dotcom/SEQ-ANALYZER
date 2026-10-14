#!/usr/bin/env python3
"""
build_database.py  -  builds reference.db from UniProt (run this ONCE, needs internet)

UniProt/Swiss-Prot is a library of expert-checked protein records. For each protein we save:
  - its sequence (the "correct copy" your sequence is compared against)
  - the diseases linked to it
  - its known natural variants (mutations), with the disease each one is linked to

Usage
  python build_database.py                      build with the starter list of genes
  python build_database.py --add P68871 P01308  add specific UniProt accessions
  python build_database.py --search "gene:BRCA1 AND organism_id:9606 AND reviewed:true"
  python build_database.py --file my_accessions.txt   (one accession per line)
  python build_database.py --check              show what is in the database now
  python build_database.py --fresh              delete the old database and start over

Uses only Python's built-in libraries. No extra install needed.
"""

import argparse
import datetime
import json
import os
import re
import sqlite3
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

DB_PATH = os.environ.get("SEQ_DB") or os.path.join(os.path.dirname(os.path.abspath(__file__)), "reference.db")
BASE = "https://rest.uniprot.org/uniprotkb"

# A starter list of well-studied human genes linked to disease. Each accession is only a
# pointer: the name, sequence and variants all come from UniProt itself when you run this.
STARTER = [
    "P68871",  # HBB    hemoglobin beta (sickle cell, beta-thalassemia)
    "P69905",  # HBA1   hemoglobin alpha
    "P01308",  # INS    insulin
    "P04637",  # TP53   tumor suppressor
    "P38398",  # BRCA1  breast/ovarian cancer
    "P51587",  # BRCA2
    "P13569",  # CFTR   cystic fibrosis
    "P01116",  # KRAS
    "P15056",  # BRAF
    "P00533",  # EGFR
    "P60484",  # PTEN
    "P04626",  # ERBB2  (HER2)
    "P42336",  # PIK3CA
    "P00519",  # ABL1
    "P10721",  # KIT
    "P07949",  # RET
    "P21359",  # NF1
    "P40692",  # MLH1   Lynch syndrome
    "P43246",  # MSH2   Lynch syndrome
    "P42858",  # HTT    Huntington disease
    "P05067",  # APP    Alzheimer disease
    "P37840",  # SNCA   Parkinson disease
    "P02649",  # APOE
    "P00439",  # PAH    phenylketonuria
    "P35498",  # SCN1A  epilepsy
    "P02768",  # ALB    albumin
    "P42212",  # GFP    green fluorescent protein (a non-human control)
]

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);
CREATE TABLE IF NOT EXISTS proteins (
    accession TEXT PRIMARY KEY, gene TEXT, name TEXT, organism TEXT, length INTEGER, sequence TEXT);
CREATE TABLE IF NOT EXISTS protein_diseases (
    accession TEXT, acronym TEXT, name TEXT, mim TEXT);
CREATE TABLE IF NOT EXISTS variants (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    accession TEXT, position INTEGER, ref TEXT, alt TEXT,
    description TEXT, diseases TEXT, dbsnp TEXT, feature_id TEXT);
CREATE INDEX IF NOT EXISTS idx_variants ON variants (accession, position);
CREATE INDEX IF NOT EXISTS idx_pd ON protein_diseases (accession);
"""

AMINO = set("ACDEFGHIKLMNPQRSTVWY")


# ---------------------------------------------------------------
# Talking to UniProt
# ---------------------------------------------------------------
def fetch_json(url, tries=3):
    req = urllib.request.Request(url, headers={
        "Accept": "application/json", "User-Agent": "sequence-analyzer-student-project"})
    last = None
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return json.load(resp), resp.headers.get("X-UniProt-Release", "")
        except urllib.error.HTTPError as e:
            if e.code in (400, 404, 410):
                return None, ""
            last = e
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = e
        time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"Could not reach UniProt ({last}). Check your internet connection. "
                       f"If you see a certificate error, update Python or try another network.")


# ---------------------------------------------------------------
# Reading one UniProt record (written defensively: missing parts are skipped, never crash)
# ---------------------------------------------------------------
def _get(d, *path, default=""):
    for key in path:
        if isinstance(d, dict):
            d = d.get(key)
        elif isinstance(d, list) and isinstance(key, int) and 0 <= key < len(d):
            d = d[key]
        else:
            return default
        if d is None:
            return default
    return d


def parse_entry(entry):
    """Turn one UniProt JSON record into plain Python data."""
    acc = _get(entry, "primaryAccession")
    if not acc:
        return None

    sequence = str(_get(entry, "sequence", "value")).upper()
    sequence = "".join(c if c in AMINO else "X" for c in sequence)
    if not sequence:
        return None

    gene = _get(entry, "genes", 0, "geneName", "value") or _get(entry, "uniProtkbId")
    name = (_get(entry, "proteinDescription", "recommendedName", "fullName", "value")
            or _get(entry, "proteinDescription", "submissionNames", 0, "fullName", "value")
            or _get(entry, "uniProtkbId"))
    organism = _get(entry, "organism", "scientificName")

    diseases = []
    for comment in _get(entry, "comments", default=[]) or []:
        if _get(comment, "commentType") != "DISEASE":
            continue
        disease = _get(comment, "disease", default={}) or {}
        mim = ""
        xref = _get(disease, "diseaseCrossReference", default={}) or {}
        if _get(xref, "database") == "MIM":
            mim = str(_get(xref, "id"))
        if _get(disease, "diseaseId"):
            diseases.append({"acronym": _get(disease, "acronym"),
                             "name": _get(disease, "diseaseId"), "mim": mim})

    variants = []
    for feat in _get(entry, "features", default=[]) or []:
        if _get(feat, "type") != "Natural variant":
            continue
        start = _get(feat, "location", "start", "value", default=None)
        if not isinstance(start, int):
            continue
        description = str(_get(feat, "description"))
        found = []
        for d in diseases:
            acr = d["acronym"]
            if acr and re.search(r"(?<![A-Za-z0-9])" + re.escape(acr) + r"(?![A-Za-z0-9])", description):
                found.append(d["name"])
        dbsnp = ""
        for xr in _get(feat, "featureCrossReferences", default=[]) or []:
            if _get(xr, "database") == "dbSNP":
                dbsnp = str(_get(xr, "id"))
                break
        original = str(_get(feat, "alternativeSequence", "originalSequence")).upper()
        alts = _get(feat, "alternativeSequence", "alternativeSequences", default=[]) or []
        if not original:
            continue
        for alt in (alts or [""]):          # no alternative = the residues are missing (deletion)
            variants.append({
                "position": start, "ref": original, "alt": str(alt).upper(),
                "description": description, "diseases": "|".join(found),
                "dbsnp": dbsnp, "feature_id": str(_get(feat, "featureId"))})

    return {"accession": acc, "gene": gene, "name": name, "organism": organism,
            "sequence": sequence, "diseases": diseases, "variants": variants}


def store_entry(conn, parsed):
    acc = parsed["accession"]
    conn.execute("DELETE FROM proteins WHERE accession=?", (acc,))
    conn.execute("DELETE FROM protein_diseases WHERE accession=?", (acc,))
    conn.execute("DELETE FROM variants WHERE accession=?", (acc,))
    conn.execute("INSERT INTO proteins VALUES (?,?,?,?,?,?)",
                 (acc, parsed["gene"], parsed["name"], parsed["organism"],
                  len(parsed["sequence"]), parsed["sequence"]))
    for d in parsed["diseases"]:
        conn.execute("INSERT INTO protein_diseases VALUES (?,?,?,?)",
                     (acc, d["acronym"], d["name"], d["mim"]))
    for v in parsed["variants"]:
        conn.execute("INSERT INTO variants (accession, position, ref, alt, description, diseases, dbsnp, feature_id) "
                     "VALUES (?,?,?,?,?,?,?,?)",
                     (acc, v["position"], v["ref"], v["alt"], v["description"],
                      v["diseases"], v["dbsnp"], v["feature_id"]))


def set_meta(conn, release):
    rows = {"built": datetime.datetime.now().strftime("%Y-%m-%d %H:%M"),
            "source": "UniProtKB (rest.uniprot.org)"}
    if release:
        rows["uniprot_release"] = release
    for k, v in rows.items():
        conn.execute("INSERT OR REPLACE INTO meta VALUES (?,?)", (k, v))


def summary_line(p):
    dv = sum(1 for v in p["variants"] if v["diseases"])
    return (f"{p['accession']:<8} {p['gene'][:10]:<10} {p['name'][:42]:<42} "
            f"{len(p['sequence']):>5} aa  {len(p['variants']):>4} variants ({dv} with disease)")


# ---------------------------------------------------------------
# Commands
# ---------------------------------------------------------------
def add_accessions(conn, accessions):
    release, ok, failed = "", 0, []
    for acc in accessions:
        data, rel = fetch_json(f"{BASE}/{urllib.parse.quote(acc)}.json")
        release = rel or release
        parsed = parse_entry(data) if data else None
        if not parsed:
            print(f"  skipped {acc}: not found or unreadable")
            failed.append(acc)
            continue
        store_entry(conn, parsed)
        conn.commit()
        print("  " + summary_line(parsed))
        ok += 1
    set_meta(conn, release)
    conn.commit()
    return ok, failed


def add_search(conn, query, size):
    url = f"{BASE}/search?" + urllib.parse.urlencode({"query": query, "format": "json", "size": size})
    data, release = fetch_json(url)
    ok = 0
    for entry in (data or {}).get("results", []):
        parsed = parse_entry(entry)
        if parsed:
            store_entry(conn, parsed)
            print("  " + summary_line(parsed))
            ok += 1
    set_meta(conn, release)
    conn.commit()
    return ok


def check(path):
    if not os.path.exists(path):
        print("No database yet. Run: python build_database.py")
        return
    conn = sqlite3.connect(path)
    meta = dict(conn.execute("SELECT key, value FROM meta").fetchall())
    n = conn.execute("SELECT COUNT(*) FROM proteins").fetchone()[0]
    v = conn.execute("SELECT COUNT(*) FROM variants").fetchone()[0]
    dv = conn.execute("SELECT COUNT(*) FROM variants WHERE diseases != ''").fetchone()[0]
    print(f"Database: {path}")
    print(f"  built: {meta.get('built', '?')}   UniProt release: {meta.get('uniprot_release', '?')}")
    print(f"  {n} proteins, {v} variants ({dv} linked to a disease)")
    print("  Proteins:")
    for row in conn.execute("SELECT accession, gene, name, length FROM proteins ORDER BY gene"):
        print(f"    {row[0]:<8} {row[1]:<10} {row[2][:50]:<50} {row[3]} aa")
    sample = conn.execute("SELECT accession, position, ref, alt, diseases, description FROM variants "
                          "WHERE diseases != '' LIMIT 3").fetchall()
    if sample:
        print("  Sample disease variants (so you can check they look right):")
        for s in sample:
            print(f"    {s[0]} {s[2]}{s[1]}{s[3] or 'del'}  ->  {s[4]}")
    conn.close()


def main():
    ap = argparse.ArgumentParser(description="Build the reference database from UniProt.")
    ap.add_argument("--add", nargs="+", metavar="ACC", help="UniProt accessions to add")
    ap.add_argument("--file", help="text file with one accession per line")
    ap.add_argument("--search", help="a UniProt search query")
    ap.add_argument("--size", type=int, default=50, help="max results for --search (default 50)")
    ap.add_argument("--check", action="store_true", help="show what is in the database")
    ap.add_argument("--fresh", action="store_true", help="delete the old database first")
    ap.add_argument("--db", default=DB_PATH, help="database file (default: reference.db)")
    args = ap.parse_args()

    if args.check:
        check(args.db)
        return

    if args.fresh and os.path.exists(args.db):
        os.remove(args.db)

    accessions = []
    if args.add:
        accessions += args.add
    if args.file:
        with open(args.file, encoding="utf8") as fh:
            accessions += [ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")]
    if not accessions and not args.search:
        accessions = STARTER

    conn = sqlite3.connect(args.db)
    conn.executescript(SCHEMA)
    print(f"Building {args.db} ...")
    try:
        if accessions:
            ok, failed = add_accessions(conn, accessions)
            print(f"\nAdded {ok} protein(s)." + (f" Skipped: {', '.join(failed)}" if failed else ""))
        if args.search:
            ok = add_search(conn, args.search, args.size)
            print(f"\nAdded {ok} protein(s) from the search.")
    except RuntimeError as err:
        print(f"\nError: {err}")
        sys.exit(1)
    finally:
        conn.close()
    print("Done. Run  python build_database.py --check  to see what is inside.")


if __name__ == "__main__":
    main()
