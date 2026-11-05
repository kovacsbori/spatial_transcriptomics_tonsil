import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from tonsil_ui import (REGION_COLORS, add_value_trace, gene_values, how_it_works, in_crop, question, region_crop,
                       results, tissue_fig, what_you_see)

R = results()
adata = R["adata"]
regions = adata.obs["region"].astype(str).values
nh = R["nhood"]
moran = R["moran"]

question("Which tissue parts are neighbours, and which keep apart?")

# ------------------------------------------------------------------------------------- part 1: neighbours
st.subheader("Part 1 · Counting neighbours")
what_you_see(
    "**Left:** every spot has up to 6 direct neighbours (black lines). **Right:** for each pair of tissue parts, "
    "do they touch **more** (red) or **less** (blue) than if the spots were shuffled randomly?"
)
c1, c2 = st.columns([1, 1.5])
with c1:
    zoom = st.segmented_control("Zoom in on", ["Germinal center", "T-cell zone", "Squamous epithelium"],
                                default="Germinal center", required=True)
    crop = region_crop(zoom, half_width=90)
    idx = np.where(in_crop(crop))[0]
    A = adata.obsp["spatial_connectivities"][idx][:, idx].tocoo()
    gfig = tissue_fig(height=420, crop=crop, force_image=True)
    xs, ys = [], []
    for a, b in zip(A.row, A.col):
        if a < b:
            xs += [R["xy"][idx[a], 0], R["xy"][idx[b], 0], None]
            ys += [R["xy"][idx[a], 1], R["xy"][idx[b], 1], None]
    gfig.add_trace(go.Scatter(x=xs, y=ys, mode="lines", line=dict(color="black", width=1.5), hoverinfo="skip",
                              showlegend=False))
    for reg, color in REGION_COLORS.items():
        mm = regions[idx] == reg
        if mm.any():
            gfig.add_trace(go.Scatter(x=R["xy"][idx[mm], 0], y=R["xy"][idx[mm], 1], mode="markers", name=reg,
                                      marker=dict(size=11, color=color, line=dict(width=1, color="white"))))
    gfig.update_layout(legend=dict(orientation="h", y=-0.02))
    st.plotly_chart(gfig, alt=f"Zoomed view of neighbouring spots around a {zoom.lower()}")
with c2:
    hfig = px.imshow(nh.round(0), text_auto=True, color_continuous_scale="RdBu_r", zmin=-40, zmax=40, aspect="auto")
    hfig.update_layout(height=470, margin=dict(l=0, r=0, t=10, b=0),
                       coloraxis_colorbar=dict(title="", tickvals=[-35, 0, 35],
                                               ticktext=["keep apart", "random", "touch"]))
    hfig.update_xaxes(tickangle=-35)
    st.plotly_chart(hfig, alt="Heatmap showing which tissue parts touch more or less than chance")
    st.caption("Ignore the diagonal (a region always touches itself). The interesting squares are off the diagonal.")

# ------------------------------------------------------------------------------------- part 2: patterned genes
st.space("large")
st.subheader("Part 2 · Genes that form patterns")
what_you_see(
    "We asked the computer, without any labels: *which genes form clear patches?* Here are some of the top hits."
)
genes = [g for g in ["CCL21", "CCL19", "FDCSP", "RGS13", "S100A8", "SPRR3"] if g in moran.index]
gene = st.pills("Gene", genes, default=genes[0], required=True)
c1, c2 = st.columns([1.5, 1])
with c1:
    mfig = tissue_fig(height=460)
    add_value_trace(mfig, gene_values(gene), gene)
    st.plotly_chart(mfig, alt=f"Tissue map of {gene}")
with c2:
    with st.container(border=True):
        st.markdown(
            "- *CCL21*, *CCL19*: signals that **attract T cells** and keep them in the T-cell zone\n"
            "- *FDCSP*, *RGS13*: **germinal-centre** genes\n"
            "- *S100A8*, *SPRR3*: genes of the **surface lining**"
        )
    st.metric("Pattern strength (0 = random, 1 = perfect patches)", f"{moran.loc[gene, 'I']:.2f}", border=True)

how_it_works([
    "Spots sit on a honeycomb grid, so each has up to 6 neighbours. These links are the 'spatial graph'.",
    "Neighbourhood enrichment: count how often two tissue parts touch, compare with 1,000 random shuffles of the "
    "labels. The number is a z-score: positive = more contact than chance.",
    "Pattern strength is Moran's I: does a spot's value resemble its neighbours' values?",
])
