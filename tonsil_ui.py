"""Shared helpers for the tonsil walkthrough app: cached results, tissue plots and page layout blocks.

Every page in app_pages/ imports from here, so all slides look and behave the same.
"""
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st
from scipy import sparse

from tonsil_analysis import run

# Fixed colour per named tissue region, used on every slide so a region always has the same colour.
REGION_COLORS = {
    "B follicle (mantle)": "#4C78A8",
    "Germinal center": "#E45756",
    "T-cell zone": "#54A24B",
    "Plasma cells": "#B279A2",
    "Squamous epithelium": "#F2B701",
    "Stroma / vessels": "#9D755D",
}
# Palette for numbered clusters (the clustering slide can produce up to ~18 clusters).
CLUSTER_PALETTE = px.colors.qualitative.Dark24


# ------------------------------------------------------------------------------------------- data access
@st.cache_resource(show_spinner="Running the full analysis once (~25 s)…")
def results():
    """Run the pipeline once per server and share the result with every page and every viewer.

    cache_resource keeps the returned objects in memory (no copying), which is what we want for the large
    AnnData object. The first page load takes ~25 s; afterwards every call returns instantly.
    """
    R = run()
    # Sparse matrices are stored row-wise (CSR) = fast per spot. Plots need one gene across all spots
    # (a column), so we keep column-wise (CSC) copies for fast gene lookup.
    R["X_csc"] = R["adata"].X.tocsc()                       # log-normalised expression
    R["counts_csc"] = R["adata"].layers["counts"].tocsc()   # raw counts
    return R


def gene_values(g, layer="norm"):
    """Expression of gene `g` in every spot, as a 1-D array (layer: "norm" = log-normalised, else raw counts)."""
    R = results()
    M = R["X_csc"] if layer == "norm" else R["counts_csc"]
    col = M[:, R["adata"].var_names.get_loc(g)]
    return np.asarray(col.toarray() if sparse.issparse(col) else col).ravel()


@st.cache_data(max_entries=20, show_spinner="Re-clustering…")
def leiden_at(resolution):
    """Leiden clusters at another resolution, for the 'level of detail' slider (labels as strings).

    Reuses the stored expression neighbour graph, so only the (fast) community detection is redone.
    Cached per resolution, so moving the slider back and forth is instant.
    """
    import anndata as ad
    import pandas as pd
    import scanpy as sc

    from tonsil_analysis import SEED

    R = results()
    # Leiden writes its result into an AnnData; use an empty one with the right spots so the shared
    # (cached) object is never modified.
    tmp = ad.AnnData(obs=pd.DataFrame(index=R["adata"].obs_names))
    sc.tl.leiden(tmp, resolution=resolution, adjacency=R["adata"].obsp["connectivities"], flavor="igraph",
                 n_iterations=2, directed=False, random_state=SEED)
    return tmp.obs["leiden"].astype(str).values


def display_opts():
    """Current sidebar display settings: (spot size, spot opacity, show H&E image)."""
    s = st.session_state
    return s.get("spot_size", 6), s.get("opacity", 0.85), s.get("show_he", True)


# ------------------------------------------------------------------------------------------- figures
def tissue_fig(height=600, crop=None, force_image=False):
    """Empty Plotly figure set up as a tissue map: H&E image in the background, axes in hires-image pixels.

    crop = (x0, x1, y0, y1) shows only part of the image (used for zoomed views).
    force_image=True shows the H&E even if the sidebar toggle is off (for slides that are about the image).
    Spots are added on top with the add_*_trace helpers below.
    """
    R = results()
    _, _, show_he = display_opts()
    if show_he or force_image:
        img = R["image"]
        x0 = y0 = 0
        if crop:
            x0, x1, y0, y1 = (int(v) for v in crop)
            img = img[y0:y1, x0:x1]
        # binary_string + jpg: send the image as one compressed JPEG instead of millions of pixel values.
        # x/y offsets keep the coordinates of a cropped image equal to the full image's coordinates.
        fig = px.imshow(img, binary_string=True, binary_format="jpg", x=np.arange(img.shape[1]) + x0,
                        y=np.arange(img.shape[0]) + y0)
        fig.data[0].hoverinfo = "skip"          # hovering should show spot info, not pixel colours
    else:
        fig = go.Figure()
        # Image convention: y grows downwards; keep x and y on the same scale so the tissue isn't stretched.
        fig.update_yaxes(autorange="reversed", scaleanchor="x")
    fig.update_xaxes(visible=False)
    fig.update_yaxes(visible=False)
    # Transparent background so the plot blends into the dark theme.
    fig.update_layout(height=height, margin=dict(l=0, r=0, t=30, b=0), legend=dict(itemsizing="constant"),
                      plot_bgcolor="rgba(0,0,0,0)")
    return fig


def add_category_traces(fig, x, y, labels, colors, highlight=None, size=None, name_prefix=""):
    """Add one scatter trace per category (e.g. region or cluster), so each gets its own legend entry.

    highlight: optional boolean mask; spots outside it are drawn dim grey (used for lasso selections).
    """
    spot_size, opacity, _ = display_opts()
    size = size or spot_size
    for cat, color in colors.items():
        m = labels == cat
        if highlight is not None:
            # Unselected spots of this category: grey, no legend, no hover
            lo = m & ~highlight
            fig.add_trace(go.Scattergl(x=x[lo], y=y[lo], mode="markers", showlegend=False, hoverinfo="skip",
                                       marker=dict(size=size, color="#475569", opacity=0.45)))
            m = m & highlight
        # legendgroup links traces with the same category across panels (clicking the legend toggles all)
        fig.add_trace(go.Scattergl(x=x[m], y=y[m], mode="markers", name=f"{name_prefix}{cat}", legendgroup=str(cat),
                                   marker=dict(size=size, color=color, opacity=opacity),
                                   hovertemplate=f"{name_prefix}{cat}<extra></extra>"))
    return fig


def add_region_traces(fig, highlight=None):
    """Spots coloured by their named tissue region (the most common view in the app)."""
    R = results()
    regions = R["adata"].obs["region"].astype(str).values
    return add_category_traces(fig, R["xy"][:, 0], R["xy"][:, 1], regions, REGION_COLORS, highlight)


def add_value_trace(fig, values, label, cmap="Viridis", x=None, y=None, hover_extra=None):
    """Spots coloured by a continuous value (gene expression, score, protein …) with a colour bar.

    The colour range is clipped to the 1st–99th percentile so a few extreme spots don't wash out the map.
    x/y default to the filtered spots' positions; pass them explicitly for other spot sets (e.g. QC_all).
    """
    R = results()
    spot_size, opacity, _ = display_opts()
    x = R["xy"][:, 0] if x is None else x
    y = R["xy"][:, 1] if y is None else y
    hover = hover_extra if hover_extra is not None else R["adata"].obs["region"].astype(str).values
    fig.add_trace(go.Scattergl(
        x=x, y=y, mode="markers", showlegend=False,
        marker=dict(size=spot_size, color=values, colorscale=cmap, opacity=opacity,
                    # tiny epsilon keeps cmax > cmin even for a gene that is constant across spots
                    cmin=np.percentile(values, 1), cmax=max(np.percentile(values, 99), np.percentile(values, 1) + 1e-6),
                    colorbar=dict(title=label, thickness=12)),
        customdata=np.stack([hover, np.round(values, 2)], axis=1),
        hovertemplate="%{customdata[0]}<br>" + label + ": %{customdata[1]}<extra></extra>"))
    return fig


def region_crop(region, half_width=110):
    """Crop window (x0, x1, y0, y1) centred on the spot with the highest score for a region's signature.

    Used for the zoomed views, so 'Zoom in on: Germinal center' always lands on a typical germinal centre.
    """
    R = results()
    obs = R["adata"].obs
    i = int(np.argmax(obs[region].values))       # the signature score column has the region's name
    cx, cy = R["xy"][i]
    h, w = R["image"].shape[:2]
    # Keep the window inside the image
    x0, y0 = max(0, cx - half_width), max(0, cy - half_width)
    x1, y1 = min(w, x0 + 2 * half_width), min(h, y0 + 2 * half_width)
    return x0, x1, y0, y1


def in_crop(crop):
    """Boolean mask of the spots that fall inside a crop window."""
    R = results()
    x0, x1, y0, y1 = crop
    x, y = R["xy"][:, 0], R["xy"][:, 1]
    return (x >= x0) & (x < x1) & (y >= y0) & (y < y1)


# ------------------------------------------------------------------------------------------- page layout
# Small building blocks so every slide has the same structure: question -> explanation -> figure ->
# optional technical detail.
def question(text):
    """The slide's guiding question, as a heading in the accent colour."""
    st.markdown(f"### :primary[{text}]")


def what_you_see(text):
    """One or two plain sentences telling the audience what the figure shows / what to do."""
    st.markdown(text)


def how_it_works(steps):
    """Optional technical detail, collapsed by default so the page stays simple."""
    with st.expander("How it works (optional detail)", icon=":material/settings:"):
        st.markdown("\n".join(f"{i}. {s}" for i, s in enumerate(steps, 1)))


def references(*keys):
    """Short list of key references (used on the final slide). Keys refer to tonsil_refs.REFS."""
    from tonsil_refs import REFS
    with st.container(border=True):
        st.markdown("**:material/menu_book: Key references**")
        st.caption("  \n".join(f"**{REFS[k][0]}**: {REFS[k][1]}. [{REFS[k][2].replace('https://', '')}]({REFS[k][2]})"
                                for k in keys))
