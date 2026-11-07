import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from tonsil_ui import (REGION_COLORS, add_value_trace, gene_values, how_it_works, in_crop, question, region_crop,
                       results, tissue_fig)

R = results()
adata = R["adata"]
regions = adata.obs["region"].astype(str).values
nh = R["nhood"]
# Short region names for the narrow controls and heatmap axes
SHORT = {"B follicle (mantle)": "Mantle", "Germinal center": "Germ. center", "T-cell zone": "T-cell zone",
         "Plasma cells": "Plasma", "Squamous epithelium": "Epithelium", "Stroma / vessels": "Stroma"}
moran = R["moran"]

question("Which tissue parts are neighbours, and which keep apart?")

c1, c2, c3 = st.columns([1, 1, 1], gap="medium")

# ------------------------------------------------------------------------------------- part 1: neighbours
with c1:
    st.markdown("**1 · Each spot has up to 6 neighbours**")
    zoom = st.segmented_control("Zoom in on", ["Germinal center", "T-cell zone", "Squamous epithelium"],
                                default="Germinal center", required=True, label_visibility="collapsed",
                                format_func=lambda r: SHORT[r])
    crop = region_crop(zoom, half_width=90)
    idx = np.where(in_crop(crop))[0]
    A = adata.obsp["spatial_connectivities"][idx][:, idx].tocoo()
    gfig = tissue_fig(height=280, crop=crop, force_image=True)
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
                                      showlegend=False, marker=dict(size=10, color=color,
                                                                    line=dict(width=1, color="white"))))
    st.plotly_chart(gfig, alt=f"Zoomed view of neighbouring spots around a {zoom.lower()}")
    how_it_works([
        "Spatial graph: hexagonal Visium grid, 6 direct neighbours per spot.",
        "Neighbourhood enrichment: contacts between two regions vs. 1,000 label shuffles → z-score "
        "(+ touch more than chance, − avoid).",
        "Pattern strength = Moran's I: similarity of a gene's value to its neighbours' values (0 random, 1 patches).",
    ])

with c2:
    st.markdown("**2 · Who touches whom?**  \n:gray[Red = touch more than chance, blue = keep apart. "
                "Ignore the diagonal.]")
    # Short region names on the axes, so the labels don't overlap in the narrow column
    short = nh.rename(index=SHORT, columns=SHORT).round(0)
    hfig = px.imshow(short, text_auto=True, color_continuous_scale="RdBu_r", zmin=-40, zmax=40, aspect="auto",
                     labels=dict(x="", y="", color="z"))
    hfig.update_traces(textfont_size=13)
    # No colour bar: the caption above explains the colours, and the heatmap gets the full column width
    hfig.update_layout(height=400, margin=dict(l=0, r=0, t=10, b=0), coloraxis_showscale=False)
    hfig.update_xaxes(tickangle=-40, side="bottom")
    st.plotly_chart(hfig, alt="Heatmap showing which tissue parts touch more or less than chance")

# ------------------------------------------------------------------------------------- part 3: patterned genes
with c3:
    st.markdown("**3 · Genes that form patches**  \n:gray[Top hits found without any labels.]")
    genes = [g for g in ["CCL21", "CCL19", "FDCSP", "RGS13", "S100A8", "SPRR3"] if g in moran.index]
    gene = st.pills("Gene", genes, default=genes[0], required=True, label_visibility="collapsed")
    mfig = tissue_fig(height=400, title=f"{gene} · pattern strength {moran.loc[gene, 'I']:.2f}")
    add_value_trace(mfig, gene_values(gene), gene)
    st.plotly_chart(mfig, alt=f"Tissue map of {gene}")
    st.caption("*CCL21/CCL19*: attract T cells to the T-cell zone · *FDCSP/RGS13*: germinal centre · "
               "*S100A8/SPRR3*: surface lining")
