import numpy as np
import pandas as pd
import plotly.express as px
import streamlit as st

from tonsil_ui import add_value_trace, how_it_works, question, results, tissue_fig, what_you_see

R = results()
adata = R["adata"]
counts = R["counts_csc"].tocsr()
obs = adata.obs

@st.dialog("The data as a table", width="large")
def table_dialog(i):
    """Large pop-up with example rows of the spot × gene table (the clicked spot first, if any)."""
    genes = [g for g in ["MS4A1", "CD3E", "RGS13", "KRT13", "JCHAIN", "CCL21"] if g in adata.var_names]
    others = np.random.default_rng(0).choice(adata.n_obs, 7 if i is not None else 8, replace=False)
    sample = others if i is None else np.r_[i, others]
    cols = [adata.var_names.get_loc(g) for g in genes]
    table = pd.DataFrame(counts[sample][:, cols].toarray().astype(int), columns=genes,
                         index=[f"spot {k + 1}" for k in range(len(sample))])
    table.insert(0, "y position", R["xy"][sample, 1].round(0).astype(int))
    table.insert(0, "x position", R["xy"][sample, 0].round(0).astype(int))
    st.dataframe(table, alt="Example rows of the spot-by-gene table with positions")
    st.caption("Rows = spots, columns = genes. The position columns make it *spatial*.")


question("What do we actually measure in one spot?")
what_you_see("**Click any spot** on the tissue to see its genes. With nothing selected you see the whole tissue; "
             "double-click the map to go back.")

c1, c2, c3 = st.columns([1.15, 1, 0.8], gap="medium")
with c1:
    fig = tissue_fig(height=470)
    add_value_trace(fig, obs["total_counts"].values.astype(float), "RNA amount")
    event = st.plotly_chart(fig, on_select="rerun", selection_mode="points", key="spot_click",
                            alt="Tissue map coloured by RNA amount per spot; click a spot to inspect it")
    pts = event.selection.get("points", []) if event else []
    i = pts[-1]["point_index"] if pts else None          # None = nothing selected -> whole tissue

with c2:
    if i is None:
        # Whole tissue: sum over all spots
        gene_totals = np.asarray(counts.sum(axis=0)).ravel()
        top = pd.Series(gene_totals, index=adata.var_names).sort_values(ascending=False).head(10)
        n_molecules, n_genes = int(gene_totals.sum()), int((gene_totals > 0).sum())
        title = f"Top 10 genes · all {adata.n_obs:,} spots"
    else:
        row = counts[i]
        top = pd.Series(row.data, index=adata.var_names[row.indices]).sort_values(ascending=False).head(10)
        n_molecules, n_genes = int(obs["total_counts"].iloc[i]), int(obs["n_genes_by_counts"].iloc[i])
        title = "Top 10 genes in this spot"
    with st.container(horizontal=True):
        st.metric("RNA molecules", f"{n_molecules:,}", border=True)
        st.metric("Different genes", f"{n_genes:,}", border=True)
    bfig = px.bar(top[::-1], orientation="h", labels={"value": "molecules counted", "index": ""})
    bfig.update_layout(height=330, margin=dict(l=0, r=0, t=30, b=0), showlegend=False, title=title)
    st.plotly_chart(bfig, alt="Bar chart of the most abundant genes in the selected spot")

with c3:
    how_it_works([
        "Space Ranger (10x) assigns each read to a gene and a spot barcode.",
        "Unique molecular identifiers (UMIs) collapse PCR copies, so each molecule counts once.",
        "Only barcodes under the tissue are kept (the 'filtered' matrix).",
        f"Result: {adata.n_obs:,} spots × {adata.n_vars:,} genes, plus x/y image coordinates per spot.",
    ])
    if st.button("See the data as a table", icon=":material/table:", width="stretch"):
        table_dialog(i)
