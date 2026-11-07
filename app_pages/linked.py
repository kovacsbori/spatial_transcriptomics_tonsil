import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import streamlit as st

from tonsil_ui import REGION_COLORS, add_region_traces, how_it_works, question, results, tissue_fig, what_you_see

R = results()
adata = R["adata"]
regions = adata.obs["region"].astype(str).values
umap = adata.obsm["X_umap"]

question("Are spots with similar genes also close together in the tissue?")
what_you_see(
    "**Draw a loop** around some dots on the similarity map (left). The same spots light up on the tissue (middle). "
    "Double-click the left plot to clear."
)

ufig = go.Figure()
trace_idx = []   # maps trace number -> spot indices, to translate a selection back to spots
for reg, color in REGION_COLORS.items():
    idx = np.where(regions == reg)[0]
    trace_idx.append(idx)
    ufig.add_trace(go.Scattergl(x=umap[idx, 0], y=umap[idx, 1], mode="markers", name=reg,
                                marker=dict(size=4, color=color), hovertemplate=f"{reg}<extra></extra>"))
ufig.update_layout(height=450, margin=dict(l=0, r=0, t=30, b=0), dragmode="lasso", plot_bgcolor="rgba(0,0,0,0)",
                   legend=dict(orientation="h", y=-0.03, itemsizing="constant"),
                   xaxis=dict(visible=False), yaxis=dict(visible=False), title="Similarity map")

c1, c2, c3 = st.columns([1, 1, 0.8], gap="medium")
with c1:
    event = st.plotly_chart(ufig, on_select="rerun", selection_mode=("lasso", "box"), key="umap_sel",
                            alt="Similarity map of spots coloured by region; draw a loop to select")
selected = np.zeros(adata.n_obs, dtype=bool)
for p in (event.selection.get("points", []) if event else []):
    selected[trace_idx[p["curve_number"]][p["point_index"]]] = True

with c2:
    tfig = tissue_fig(height=450, title="The tissue")
    add_region_traces(tfig, highlight=selected if selected.any() else None)
    tfig.update_layout(showlegend=False)
    st.plotly_chart(tfig, alt="Tissue map highlighting the selected spots")

with c3:
    if selected.any():
        X = adata.X.tocsr()
        m_in = np.asarray(X[selected].mean(axis=0)).ravel()
        m_out = np.asarray(X[~selected].mean(axis=0)).ravel() if (~selected).any() else np.zeros_like(m_in)
        diff = pd.Series(m_in - m_out, index=adata.var_names).sort_values(ascending=False).head(8)
        dfig = px.bar(diff[::-1], orientation="h", labels={"value": "higher than in the rest", "index": ""})
        dfig.update_layout(height=260, margin=dict(l=0, r=0, t=30, b=0), showlegend=False,
                           title=f"What makes these {selected.sum():,} spots special")
        dfig.update_traces(marker_color="#E45756")
        st.plotly_chart(dfig, alt="Genes most increased in the selected spots")
    else:
        st.caption("Suggestion: loop the red (germinal centre) dots first.")
    how_it_works([
        "The similarity map is built from expression only (PCA → neighbour graph → UMAP).",
        "Selected spots are mapped back to their tissue coordinates.",
        "Bars: mean log-expression in the selection minus the rest of the tissue, top 8 genes.",
    ])
