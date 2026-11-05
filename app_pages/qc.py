import plotly.graph_objects as go
import streamlit as st
from scipy.stats import spearmanr

from tonsil_analysis import QC_MIN_COUNTS
from tonsil_ui import (add_value_trace, gene_values, how_it_works, question, results, tissue_fig,
                       what_you_see)

R = results()
qc = R["qc_all"]
adata = R["adata"]

question("Can we trust every spot, and can we compare spots fairly?")

# ------------------------------------------------------------------------------------- part 1: remove bad spots
st.subheader("Part 1 · Throw out spots with too little RNA")
what_you_see("Move the slider. Spots below the line turn into **red crosses**: they would be removed.")

min_counts = st.slider("Minimum RNA molecules per spot", 0, 20000, QC_MIN_COUNTS, step=250)
fail = (qc["total_counts"] < min_counts).values

c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=500)
    # Viridis has no red in it, so the red "removed" crosses stay easy to see
    add_value_trace(fig, qc["total_counts"].values.astype(float), "RNA amount", cmap="Viridis",
                    x=qc["x"].values, y=qc["y"].values, hover_extra=qc.index.values)
    if fail.any():
        fig.add_trace(go.Scattergl(x=qc["x"][fail], y=qc["y"][fail], mode="markers", name="removed",
                                   marker=dict(symbol="x", size=11, color="#FF0000",
                                               line=dict(width=1, color="white"))))
    fig.update_layout(legend=dict(x=0, y=1))
    st.plotly_chart(fig, alt="Tissue map of RNA amount per spot with removed spots as red crosses")
with c2:
    st.metric("Spots kept", f"{(~fail).sum():,}", border=True)
    st.metric("Spots removed", f"{fail.sum():,}", border=True)
    if min_counts > 10000:
        st.warning("The removed spots are not random: whole regions disappear, because some tissue simply "
                   "contains less RNA. A cut-off this strict would delete real biology.", icon=":material/warning:")

# ------------------------------------------------------------------------------------- part 2: fair comparison
st.space("large")
st.subheader("Part 2 · Make spots comparable")
what_you_see(
    "*ACTB* is a 'housekeeping' gene: almost every cell has it at a similar level. So its map should look fairly "
    "even. Switch between **Before** and **After** adjustment."
)
scale = st.segmented_control("Values", ["Before adjustment", "After adjustment"], default="Before adjustment",
                             required=True)
gene = "ACTB"
raw, norm = gene_values(gene, "counts").astype(float), gene_values(gene, "norm")
depth = adata.obs["total_counts"].values
rho_raw, rho_norm = spearmanr(raw, depth)[0], spearmanr(norm, depth)[0]

c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=480)
    add_value_trace(fig, raw if scale == "Before adjustment" else norm, gene)
    st.plotly_chart(fig, alt=f"Tissue map of {gene} {scale.lower()}")
with c2:
    st.metric("How much the gene just follows 'total RNA in the spot'",
              f"{(rho_raw if scale == 'Before adjustment' else rho_norm):.2f}",
              help="Correlation between the gene and the total RNA per spot (1 = follows it completely, 0 = not at all).",
              border=True)
    st.caption("**Before:** the map mostly shows how much RNA each spot caught.  \n"
               "**After:** every spot is scaled to the same total, so differences are real.")

how_it_works([
    "Quality check: for each spot we count RNA molecules and different genes. We remove spots with fewer than "
    "500 molecules or 250 genes, and genes seen in fewer than 10 spots.",
    "Adjustment (normalisation): each spot's counts are scaled to 10,000 in total, then log-transformed "
    "(`log(1+x)`) so very abundant genes don't dominate.",
    "We also pick the 2,000 genes that vary most across the tissue ('highly variable genes'); only those are "
    "used in the next step.",
])
