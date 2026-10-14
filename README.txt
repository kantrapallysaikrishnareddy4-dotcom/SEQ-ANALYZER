SEQUENCE ANALYZER + DATABASE MATCHING (Python + webpage)

SET UP (one time)
1. Open this folder in VS Code and open a terminal (Terminal > New Terminal).
2. Install the libraries:
       pip install flask biopython
   (If biopython will not install, skip it. Everything still works, just slower.)
3. Build the reference database. This needs internet and takes about a minute:
       python build_database.py
   It downloads real protein records and their known disease variants from UniProt.
4. Check it worked:
       python build_database.py --check
   You should see ~27 proteins and a few sample disease variants.

RUN
5. Start the server:
       python app.py
6. Open http://127.0.0.1:5000 in your browser.
7. Stop the server with Ctrl+C in the terminal.

HOW TO READ THE "DATABASE MATCHES" CARD
- It shows which known protein your sequence matches best, how identical it is,
  and every difference from the reference (e.g. p.Glu7Val).
- Each difference is looked up in the known-variant table. A red pill means
  UniProt links that exact change to a disease.
- DNA/RNA is translated in all 6 reading frames, so genes on either strand are found.
- Positions use UniProt numbering (the first Met is residue 1). Older papers often
  skip the Met, so sickle cell is "E6V" in papers but "Glu7Val" here.
- E-value: how many matches this good you would expect by chance. Smaller = stronger.

ADD MORE GENES
       python build_database.py --add P12345 Q67890        (UniProt accessions)
       python build_database.py --file my_accessions.txt   (one accession per line)
       python build_database.py --search "gene:BRCA1 AND organism_id:9606 AND reviewed:true"
       python build_database.py --fresh                    (start over)

FILES
  app.py                 the Python server
  matcher.py             alignment, mutation calling, variant lookup
  build_database.py      downloads UniProt data into reference.db
  sequence_analyzer.py   your original analyzer
  static/                the webpage (index.html, style.css, script.js)

LIMITS (please read)
- For research and learning, not for diagnosis. A difference may be a real mutation,
  a harmless common variant, or a sequencing error. Check ClinVar for clinical meaning.
- It works at the protein level. Silent DNA changes and changes outside genes are invisible.
- Only proteins in your database can match. No match does not mean "unknown sequence".
- A frameshift shows up as a match that stops early or breaks in pieces. Only the best
  piece is reported.
- Variant data comes from UniProt's curated annotations, which is not a full ClinVar copy.
- Double-clicking index.html will not work. The page needs app.py running.
