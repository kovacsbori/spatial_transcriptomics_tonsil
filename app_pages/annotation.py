import plotly.express as px
import streamlit as st

from tonsil_analysis import SIGNATURES
from tonsil_ui import add_region_traces, add_value_trace, how_it_works, question, results, tissue_fig, what_you_see

R = results()
adata = R["adata"]
obs = adata.obs

DESCRIPTIONS = {
    "B follicle (mantle)": "Round clusters of resting B cells, the cells that make antibodies.",
    "Germinal center": "The centre of a B-cell follicle, where B cells multiply and improve their antibodies.",
    "T-cell zone": "The area between follicles where T cells meet and get activated.",
    "Plasma cells": "Mature B cells that pump out antibodies.",
    "Squamous epithelium": "The skin-like surface lining of the tonsil.",
    "Stroma / vessels": "Connective tissue and blood vessels that hold everything together.",
}

question("The groups are just numbers. Which part of the tonsil is each one?")
what_you_see(
    "**Left:** pick a tissue part; red spots are where its **textbook marker genes** are switched on. "
    "**Right:** each group named after its best-matching tissue part."
)
sig = st.segmented_control("Tissue part", list(SIGNATURES), default="Germinal center", required=True,
                           label_visibility="collapsed")

c1, c2, c3 = st.columns([1, 1, 0.75])
with c1:
    fig = tissue_fig(height=430)
    add_value_trace(fig, obs[sig].values, "marker score", cmap="RdBu_r")
    fig.update_layout(title=f"Marker score: {sig}")
    st.plotly_chart(fig, alt=f"Tissue map of the {sig} marker score")
with c2:
    fig = tissue_fig(height=430)
    add_region_traces(fig)
    fig.update_layout(title="Result: named regions", legend=dict(orientation="h", y=-0.02))
    st.plotly_chart(fig, alt="Tissue map coloured by named region")
with c3:
    with st.container(border=True):
        st.markdown(f"**{sig}**")
        st.caption(DESCRIPTIONS[sig])
        st.caption("Markers: " + ", ".join(f"*{g}*" for g in SIGNATURES[sig]))
    how_it_works([
        "Gene-set scoring ([Tirosh et al. 2016](https://doi.org/10.1126/science.aad0501); "
        "[`score_genes`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.score_genes.html)): mean "
        "expression of the marker set minus a random set of expression-matched control genes.",
        "Scores are averaged per cluster and z-scored across clusters; each cluster takes the top-scoring name.",
        "Check: each cluster's own Wilcoxon marker genes "
        "([`rank_genes_groups`](https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.rank_genes_groups.html)) "
        "must agree with its name ([more on annotation](https://www.sc-best-practices.org/cellular_structure/"
        "annotation.html)).",
    ])
    with st.popover("Scores and marker genes per group", icon=":material/table:", width="stretch"):
        hfig = px.imshow(R["cluster_z"].round(1), text_auto=True, color_continuous_scale="RdBu_r", zmin=-2, zmax=2,
                         aspect="auto", labels=dict(x="", y="group", color="score"))
        hfig.update_layout(height=340, margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(hfig, alt="Heatmap of marker score per group")
        st.dataframe(R["markers"].rename(columns=lambda c: f"group {c}"), hide_index=True,
                     alt="Top marker genes for each group")
