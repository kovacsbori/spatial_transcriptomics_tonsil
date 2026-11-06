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
    "Same tissue, two measurements: the **gene (RNA)** on the left, the **protein** (measured with antibodies) in "
    "the middle. Pick a marker and compare the two maps."
)

c1, c2, c3 = st.columns([1, 1, 0.75])
with c3:
    choice = st.segmented_control("Marker", list(PAIRS), default="PAX5 (B cells)", required=True)
    p, g = PAIRS[choice]
    rho = rp.loc[p, "rho"]
    verdict = "strong" if rho >= 0.6 else "moderate" if rho >= 0.4 else "weak"
    st.metric("Agreement (0 = none, 1 = perfect)", f"{rho:.2f} · {verdict}", border=True)
    how_it_works([
        "Antibodies carry DNA tags and are sequenced with the RNA, in the same spots.",
        "Protein counts are CLR-normalised per antibody (log counts minus that antibody's mean), removing each "
        "antibody's own background level.",
        "Agreement = Spearman correlation between RNA and protein across all spots.",
    ])
    with st.popover("Agreement for all markers", icon=":material/bar_chart:", width="stretch"):
        cfig = px.bar(R["rna_prot"].sort_values("rho"), x="rho", y="protein", orientation="h",
                      labels={"rho": "agreement (Spearman correlation)", "protein": ""})
        cfig.update_layout(height=560, margin=dict(l=0, r=0, t=10, b=0))
        st.plotly_chart(cfig, alt="Agreement between RNA and protein for every marker")

with c1:
    f1 = tissue_fig(height=470)
    add_value_trace(f1, gene_values(g), g)
    f1.update_layout(title=f"Gene (RNA): {g}")
    st.plotly_chart(f1, alt=f"Tissue map of {g} RNA")
with c2:
    f2 = tissue_fig(height=470)
    add_value_trace(f2, R["prot"][p].values, p)
    f2.update_layout(title=f"Protein: {p}")
    st.plotly_chart(f2, alt=f"Tissue map of {p} protein")
