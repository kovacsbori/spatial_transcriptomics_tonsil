import plotly.express as px
import streamlit as st

from tonsil_analysis import SIGNATURES
from tonsil_ui import (add_region_traces, add_value_trace, how_it_works, question, results, tissue_fig,
                       what_you_see)

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
    "Pick a tissue part. Each has a few **textbook marker genes**; red spots are where those genes are switched on."
)

sig = st.segmented_control("Tissue part", list(SIGNATURES), default="Germinal center", required=True)
c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=500)
    add_value_trace(fig, obs[sig].values, "marker score", cmap="RdBu_r")
    st.plotly_chart(fig, alt=f"Tissue map of the {sig} marker score")
with c2:
    with st.container(border=True):
        st.markdown(f"**{sig}**")
        st.caption(DESCRIPTIONS[sig])
        st.markdown("Marker genes: " + ", ".join(f"*{g}*" for g in SIGNATURES[sig]))
        st.caption("Red = higher than usual, blue = lower.")

st.space("large")
st.subheader("Result: each group gets the name of its best-matching tissue part")
c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=500)
    add_region_traces(fig)
    st.plotly_chart(fig, alt="Tissue map coloured by named region")
with c2:
    with st.expander("Details: score per group and each group's own top genes", icon=":material/table:"):
        hfig = px.imshow(R["cluster_z"].round(1), text_auto=True, color_continuous_scale="RdBu_r", zmin=-2, zmax=2,
                         aspect="auto", labels=dict(x="", y="group", color="score"))
        hfig.update_layout(height=360, margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(hfig, alt="Heatmap of marker score per group")
        st.dataframe(R["markers"].rename(columns=lambda c: f"group {c}"), hide_index=True,
                     alt="Top marker genes for each group")

how_it_works([
    "For each tissue part we take a list of well-known marker genes.",
    "Each spot gets a score: average activity of those genes minus a random set of similar genes (the background).",
    "Each group from the previous step gets the name of the tissue part with the highest average score.",
    "A check: the genes that make each group special (table in 'Details') should match the name.",
])
