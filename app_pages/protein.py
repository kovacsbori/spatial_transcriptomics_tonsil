import plotly.express as px
import streamlit as st

from tonsil_ui import add_value_trace, gene_values, how_it_works, question, results, tissue_fig, what_you_see

R = results()
rp = R["rna_prot"].set_index("protein")

PAIRS = {
    "PAX5 (B cells)": ("PAX5", "PAX5"),
    "CD3E (T cells)": ("CD3E", "CD3E"),
    "PanCK (surface lining)": ("PanCK", "KRT5"),
    "CD68 (macrophages)": ("CD68", "CD68"),
}

question("Genes are the recipe, proteins are the dish. Do they show the same picture?")
what_you_see(
    "Same tissue, two measurements: the **gene (RNA)** on the left, the **protein** (measured with antibodies) on "
    "the right. Pick a marker and compare the two maps."
)

choice = st.segmented_control("Marker", list(PAIRS), default="PAX5 (B cells)", required=True)
p, g = PAIRS[choice]
rho = rp.loc[p, "rho"]

c1, c2 = st.columns(2)
with c1:
    f1 = tissue_fig(height=480)
    add_value_trace(f1, gene_values(g), g)
    f1.update_layout(title=f"Gene (RNA): {g}")
    st.plotly_chart(f1, alt=f"Tissue map of {g} RNA")
with c2:
    f2 = tissue_fig(height=480)
    add_value_trace(f2, R["prot"][p].values, p, cmap="Magma")
    f2.update_layout(title=f"Protein: {p}")
    st.plotly_chart(f2, alt=f"Tissue map of {p} protein")

verdict = "strong" if rho >= 0.6 else "moderate" if rho >= 0.4 else "weak"
st.metric("Agreement between the two maps (0 = none, 1 = perfect)", f"{rho:.2f} · {verdict}", border=True)

with st.expander("Agreement for all markers", icon=":material/bar_chart:"):
    cfig = px.bar(R["rna_prot"].sort_values("rho"), x="rho", y="protein", orientation="h",
                  labels={"rho": "agreement (Spearman correlation)", "protein": ""})
    cfig.update_layout(height=620, margin=dict(l=0, r=0, t=10, b=0))
    st.plotly_chart(cfig, alt="Agreement between RNA and protein for every marker")

how_it_works([
    "The antibodies carry DNA tags, so they are read out by the same sequencing, in the same spots.",
    "Protein values are adjusted per antibody (CLR), because every antibody has its own background level.",
    "Agreement = Spearman correlation across all spots: do spots that are high for the gene also tend to be high "
    "for the protein?",
])
