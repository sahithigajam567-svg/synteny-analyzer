import io

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

import synteny as sy

st.set_page_config(page_title="Comparative Genomics Synteny Analyzer", layout="wide")
st.title("🧬 Comparative Genomics Synteny Analyzer")
st.caption("Compare gene order between two genomes and find conserved (syntenic) blocks, "
           "inversions and rearrangements. Built with Biopython.")

# ---------------- Sidebar ----------------
st.sidebar.header("1. Input")
source = st.sidebar.radio("Data source", ["Demo genomes", "Upload GenBank files"])
st.sidebar.header("2. Ortholog detection")
method = st.sidebar.radio("Method", ["Protein similarity (reciprocal best hit)",
                                     "Matching gene names"])
min_score = st.sidebar.slider("Min similarity (k-mer Jaccard)", 0.01, 0.5, 0.05, 0.01,
                              disabled=method.startswith("Matching"))
st.sidebar.header("3. Synteny blocks")
min_genes = st.sidebar.slider("Min genes per block", 2, 15, 4)
max_gap = st.sidebar.slider("Max gap between genes", 1, 20, 5)

genes_a = genes_b = None
name_a, name_b = "Genome A", "Genome B"

if source == "Demo genomes":
    genes_a, genes_b = sy.make_demo()
    st.info("Using synthetic demo genomes: B has an inversion, a translocation, "
            "a deletion and insertions relative to A. Try uploading real GenBank files in the sidebar.")
else:
    c1, c2 = st.columns(2)
    fa = c1.file_uploader("Genome A (.gb / .gbk)", type=["gb", "gbk", "gbff", "genbank"])
    fb = c2.file_uploader("Genome B (.gb / .gbk)", type=["gb", "gbk", "gbff", "genbank"])
    if fa and fb:
        try:
            genes_a = sy.parse_genbank(io.StringIO(fa.getvalue().decode("utf-8", "ignore")))
            genes_b = sy.parse_genbank(io.StringIO(fb.getvalue().decode("utf-8", "ignore")))
            name_a, name_b = fa.name, fb.name
        except Exception as e:
            st.error(f"Could not parse the files: {e}")
    else:
        st.warning("Upload two GenBank files (with CDS features and translations) to begin.")

if genes_a and genes_b:
    if not genes_a or not genes_b:
        st.error("No CDS features with translations were found.")
        st.stop()
    with st.spinner("Finding orthologs and synteny blocks..."):
        if method.startswith("Matching"):
            pairs = sy.orthologs_by_name(genes_a, genes_b)
        else:
            pairs = sy.orthologs_by_similarity(genes_a, genes_b, min_score=min_score)
        blocks = sy.find_blocks([(i, j) for i, j, _ in pairs], min_genes, max_gap)

    m = st.columns(4)
    m[0].metric(f"Genes in {name_a}", len(genes_a))
    m[1].metric(f"Genes in {name_b}", len(genes_b))
    m[2].metric("Ortholog pairs", len(pairs))
    m[3].metric("Synteny blocks", len(blocks))

    # Dot plot
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=[p[0] for p in pairs], y=[p[1] for p in pairs],
                             mode="markers", marker=dict(size=4, color="lightgray"),
                             name="Ortholog pairs"))
    for n, b in enumerate(blocks, 1):
        fwd = b["direction"] == 1
        fig.add_trace(go.Scatter(
            x=[p[0] for p in b["pairs"]], y=[p[1] for p in b["pairs"]], mode="lines+markers",
            line=dict(color="royalblue" if fwd else "crimson", width=3), marker=dict(size=5),
            name=f"Block {n} ({'forward' if fwd else 'inverted'})", showlegend=n <= 12))
    fig.update_layout(height=620, xaxis_title=f"Gene order in {name_a}",
                      yaxis_title=f"Gene order in {name_b}", title="Synteny dot plot "
                      "(blue = conserved order, red = inverted)")
    st.plotly_chart(fig, use_container_width=True)

    # Tables
    rows = [{"Block": n, f"{name_a} genes": f"{b['pairs'][0][0]}–{b['pairs'][-1][0]}",
             f"{name_b} genes": f"{min(p[1] for p in b['pairs'])}–{max(p[1] for p in b['pairs'])}",
             "Orthologs": len(b["pairs"]),
             "Orientation": "forward" if b["direction"] == 1 else "inverted"}
            for n, b in enumerate(blocks, 1)]
    st.subheader("Synteny blocks")
    bdf = pd.DataFrame(rows)
    st.dataframe(bdf, use_container_width=True)
    if rows:
        st.download_button("Download blocks (CSV)", bdf.to_csv(index=False), "synteny_blocks.csv")

    with st.expander("Ortholog pairs"):
        odf = pd.DataFrame([{"A index": i, "A gene": genes_a[i]["name"], "B index": j,
                             "B gene": genes_b[j]["name"], "Score": round(s, 3)}
                            for i, j, s in pairs])
        st.dataframe(odf, use_container_width=True)
        st.download_button("Download pairs (CSV)", odf.to_csv(index=False), "ortholog_pairs.csv")

with st.expander("How it works"):
    st.markdown("""
1. **Parse** CDS features and protein translations from GenBank files with `Bio.SeqIO`.
2. **Find orthologs** by reciprocal best hit on protein k-mer similarity (or by gene name).
3. **Chain** collinear ortholog pairs into **synteny blocks**; a run that moves backwards in
   genome B is an **inversion**. Gaps up to the chosen size are tolerated (indels).
4. **Visualise** as a dot plot and a block table.
""")
