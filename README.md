# Comparative Genomics Synteny Analyzer

Web app (Streamlit + Biopython) that compares gene order between two genomes and reports
conserved synteny blocks, inversions and rearrangements.

## Run locally (Mac)
```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```
Opens at http://localhost:8501. Works immediately with built-in demo genomes; or upload two
GenBank (.gb/.gbk) files containing CDS features with translations.

## Live demo
<paste your Streamlit Community Cloud link here>

## Files
- `app.py` – Streamlit interface
- `synteny.py` – parsing, ortholog detection, block chaining, demo data
