import plotly.graph_objects as go
import streamlit as st
from scipy.stats import spearmanr

from tonsil_analysis import QC_MIN_COUNTS
from tonsil_ui import (
    add_value_trace,
    gene_values,
    how_it_works,
    question,
    results,
    tissue_fig,
)

R = results()
qc = R["qc_all"]
adata = R["adata"]

question("Can we trust every spot, and can we compare spots fairly?")

c1, c2, c3 = st.columns([1, 1, 0.8], gap="medium")

# ------------------------------------------------------------------------------------- part 1: remove bad spots
with c1:
    st.markdown(
        "**1 · Throw out spots with too little RNA**  \n"
        ":gray[Move the slider: removed spots turn into red crosses.]"
    )
    min_counts = st.slider(
        "Minimum RNA molecules per spot",
        0,
        20000,
        QC_MIN_COUNTS,
        step=250,
        label_visibility="collapsed",
    )
    fail = (qc["total_counts"] < min_counts).values
    fig = tissue_fig(height=420)
    # Cividis has no red (or purple) in it, so the red "removed" crosses stay easy to see
    add_value_trace(
        fig,
        qc["total_counts"].values.astype(float),
        "RNA amount",
        cmap="Cividis",
        x=qc["x"].values,
        y=qc["y"].values,
        hover_extra=qc.index.values,
    )
    if fail.any():
        fig.add_trace(
            go.Scattergl(
                x=qc["x"][fail],
                y=qc["y"][fail],
                mode="markers",
                name="removed",
                marker=dict(
                    symbol="x",
                    size=11,
                    color="#FF0000",
                    line=dict(width=1, color="white"),
                ),
            )
        )
    # Legend above the map, so the "removed" key doesn't cover the tissue
    fig.update_layout(legend=dict(orientation="h", x=0, y=1, yanchor="bottom"), margin_t=26)
    st.plotly_chart(
        fig, alt="Tissue map of RNA amount per spot with removed spots as red crosses"
    )
    st.caption(
        f"Kept **{(~fail).sum():,}** · removed **{fail.sum():,}** spots (cut-off {min_counts:,})"
    )

# ------------------------------------------------------------------------------------- part 2: fair comparison
gene = "ACTB"
raw, norm = gene_values(gene, "counts").astype(float), gene_values(gene, "norm")
depth = adata.obs["total_counts"].values
rho_raw, rho_norm = spearmanr(raw, depth)[0], spearmanr(norm, depth)[0]

with c2:
    st.markdown(
        "**2 · Make spots comparable**  \n"
        ":gray[*ACTB* is in nearly every cell at a similar level, so a fair map looks even.]"
    )
    scale = st.segmented_control(
        "Values",
        ["Before adjustment", "After adjustment"],
        default="Before adjustment",
        required=True,
        label_visibility="collapsed",
    )
    fig = tissue_fig(height=420)
    add_value_trace(fig, raw if scale == "Before adjustment" else norm, gene)
    st.plotly_chart(fig, alt=f"Tissue map of {gene} {scale.lower()}")
    rho = rho_raw if scale == "Before adjustment" else rho_norm
    st.caption(
        f"Link between *ACTB* and total RNA per spot: **{rho:.2f}** (1 = follows it completely). "
        + (
            "Mostly shows how much RNA each spot caught."
            if scale == "Before adjustment"
            else "The capture effect is largely gone."
        )
    )

with c3:
    how_it_works(
        [
            f"QC: spots with < {QC_MIN_COUNTS} UMIs or < 250 genes, and genes in < 10 spots, are removed. Thresholds "
            "are mild because low RNA is often biology (fewer cells), not failure.",
            "These cut-offs are common conventions, not fixed rules; they are chosen per "
            "dataset. Here QC removes 3 spots and 55 near-empty genes.",
            "Normalisation: spots capture different total amounts of RNA, so raw counts mix biology with "
            "capture efficiency; without it, spots can't be compared fairly "
            "([more](https://www.sc-best-practices.org/preprocessing_visualization/normalization.html)).",
        ]
    )
