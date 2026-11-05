import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from tonsil_ui import add_value_trace, how_it_works, question, results, tissue_fig, what_you_see

R = results()
adata = R["adata"]
counts = R["counts_csc"].tocsr()
obs = adata.obs

question("What do we actually measure in one spot?")
what_you_see("**Click any spot** on the tissue. On the right you see the genes with the most RNA in that spot.")

c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=540)
    add_value_trace(fig, obs["total_counts"].values.astype(float), "RNA amount", cmap="Magma")
    event = st.plotly_chart(fig, on_select="rerun", selection_mode="points", key="spot_click",
                            alt="Tissue map coloured by RNA amount per spot; click a spot to inspect it")
    pts = event.selection.get("points", []) if event else []
    if pts:
        st.session_state.spot_idx = pts[-1]["point_index"]
    i = st.session_state.get("spot_idx", int(np.argmax(obs["Germinal center"].values)))

with c2:
    row = counts[i]
    top = pd.Series(row.data, index=adata.var_names[row.indices]).sort_values(ascending=False).head(10)
    with st.container(horizontal=True):
        st.metric("RNA molecules", f"{int(obs['total_counts'].iloc[i]):,}", border=True)
        st.metric("Different genes", f"{int(obs['n_genes_by_counts'].iloc[i]):,}", border=True)
    bfig = px.bar(top[::-1], orientation="h", labels={"value": "molecules counted", "index": ""})
    bfig.update_layout(height=360, margin=dict(l=0, r=0, t=30, b=0), showlegend=False,
                       title="Top 10 genes in this spot")
    st.plotly_chart(bfig, alt="Bar chart of the most abundant genes in the selected spot")

with st.expander("See the data as a table", icon=":material/table:"):
    genes = [g for g in ["MS4A1", "CD3E", "RGS13", "KRT13", "JCHAIN", "CCL21"] if g in adata.var_names]
    sample = np.r_[i, np.random.default_rng(0).choice(adata.n_obs, 7, replace=False)]
    cols = [adata.var_names.get_loc(g) for g in genes]
    table = pd.DataFrame(counts[sample][:, cols].toarray().astype(int), columns=genes,
                         index=[f"spot {k + 1}" for k in range(len(sample))])
    table.insert(0, "y position", R["xy"][sample, 1].round(0).astype(int))
    table.insert(0, "x position", R["xy"][sample, 0].round(0).astype(int))
    st.dataframe(table, alt="Example rows of the spot-by-gene table with positions")
    st.caption(f"The full table has {adata.n_obs:,} rows (spots) and {adata.n_vars:,} columns (genes). "
               "The position columns are what make it *spatial*.")

how_it_works([
    "Each sequenced read has a gene label and a spot 'postcode'. Software from 10x (Space Ranger) counts "
    "molecules per gene per spot.",
    "Copies made during sample preparation are recognised by a unique molecular tag (UMI), so each original "
    "molecule is counted once.",
    "Only spots that lie on the tissue are kept.",
])
