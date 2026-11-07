import plotly.graph_objects as go
import streamlit as st

from tonsil_analysis import LEIDEN_RESOLUTION
from tonsil_ui import (
    CLUSTER_PALETTE,
    add_category_traces,
    how_it_works,
    leiden_at,
    question,
    results,
    tissue_fig,
    what_you_see,
)

R = results()
adata = R["adata"]

question("If we sort spots only by their genes, without their position, what happens?")
what_you_see(
    "Spots with similar genes get the same colour. **Left:** a 'similarity map' where similar spots are close. "
    "**Right:** the same spots on the tissue. The computer **never saw the positions**."
)

c1, c2, c3 = st.columns([1, 1, 0.8], gap="medium")
with c3:
    res = st.select_slider(
        "Level of detail",
        options=[0.2, 0.3, 0.4, 0.6, 0.9, 1.2],
        value=LEIDEN_RESOLUTION,
        # Leiden resolution -> label with the number of groups it produces in this dataset
        format_func=lambda r: {
            0.2: "coarse · 5",
            0.3: "6",
            0.4: "7",
            0.6: "medium · 8",
            0.9: "fine · 14",
            1.2: "very fine · 18",
        }[r],
    )
    labels = leiden_at(res)
    cats = sorted(set(labels), key=int)
    colors = {c: CLUSTER_PALETTE[int(c) % len(CLUSTER_PALETTE)] for c in cats}
    st.metric("Groups found", len(cats), border=True)
    how_it_works(
        [
            "Feature selection: the 2,000 highly variable genes "
            "([seurat_v3](https://scanpy.readthedocs.io/en/stable/generated/scanpy.pp.highly_variable_genes.html)) "
            "go into clustering.",
            "UMAP draws a map where spots with similar genes sit close together; only closeness matters, not the "
            "axes or the gaps between islands ([visual guide](https://pair-code.github.io/understanding-umap/)).",
        ]
    )

umap = adata.obsm["X_umap"]
with c1:
    ufig = go.Figure()
    add_category_traces(
        ufig, umap[:, 0], umap[:, 1], labels, colors, size=4, name_prefix="group "
    )
    ufig.update_layout(
        height=450,
        margin=dict(l=0, r=0, t=30, b=0),
        title="Similarity map (UMAP)",
        xaxis=dict(visible=False),
        yaxis=dict(visible=False),
        plot_bgcolor="rgba(0,0,0,0)",
        legend=dict(itemsizing="constant"),
    )
    st.plotly_chart(ufig, alt="Similarity map of spots coloured by group")
with c2:
    tfig = tissue_fig(height=450, title="Same groups on the tissue")
    add_category_traces(
        tfig, R["xy"][:, 0], R["xy"][:, 1], labels, colors, name_prefix="group "
    )
    tfig.update_layout(showlegend=False)
    st.plotly_chart(tfig, alt="Tissue map coloured by group")
