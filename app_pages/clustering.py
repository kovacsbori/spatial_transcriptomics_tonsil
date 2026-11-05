import plotly.graph_objects as go
import streamlit as st

from tonsil_analysis import LEIDEN_RESOLUTION
from tonsil_ui import (CLUSTER_PALETTE, add_category_traces, how_it_works, leiden_at, question, results, tissue_fig, what_you_see)

R = results()
adata = R["adata"]

question("If we sort spots only by their genes, without their position, what happens?")
what_you_see(
    "The computer puts spots with similar genes into the same group (same colour). **Left:** a 'similarity map' "
    "where similar spots are close. **Right:** the same spots on the tissue. The computer **never saw the "
    "positions**."
)

res = st.select_slider("Level of detail", options=[0.2, 0.4, 0.6, 0.9, 1.2],
                       value=LEIDEN_RESOLUTION,
                       format_func=lambda r: {0.2: "coarse", 0.4: "medium-coarse", 0.6: "medium",
                                              0.9: "fine", 1.2: "very fine"}[r])
labels = leiden_at(res)
cats = sorted(set(labels), key=int)
colors = {c: CLUSTER_PALETTE[int(c) % len(CLUSTER_PALETTE)] for c in cats}
st.metric("Groups found", len(cats), border=True, width=180)

c1, c2 = st.columns(2)
umap = adata.obsm["X_umap"]
with c1:
    ufig = go.Figure()
    add_category_traces(ufig, umap[:, 0], umap[:, 1], labels, colors, size=4, name_prefix="group ")
    ufig.update_layout(height=500, margin=dict(l=0, r=0, t=30, b=0), title="Similarity map (UMAP)",
                       xaxis=dict(visible=False), yaxis=dict(visible=False), plot_bgcolor="rgba(0,0,0,0)",
                       legend=dict(itemsizing="constant"))
    st.plotly_chart(ufig, alt="Similarity map of spots coloured by group")
with c2:
    tfig = tissue_fig(height=500)
    add_category_traces(tfig, R["xy"][:, 0], R["xy"][:, 1], labels, colors, name_prefix="group ")
    tfig.update_layout(title="Same groups on the tissue", showlegend=False)
    st.plotly_chart(tfig, alt="Tissue map coloured by group")

how_it_works([
    "The 2,000 most variable genes are summarised into 30 numbers per spot (PCA), which keeps the main patterns "
    "and drops noise.",
    "Each spot is linked to the 15 spots with the most similar genes.",
    "The Leiden algorithm finds tightly linked communities = the groups ('clusters'). 'Level of detail' is "
    "its *resolution* setting.",
    "UMAP draws these links as a 2-D picture (the similarity map). Only closeness inside a group is meaningful, "
    "not the distance between far-apart groups.",
])
