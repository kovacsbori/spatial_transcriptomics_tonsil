import streamlit as st

from tonsil_ui import add_region_traces, references, results, tissue_fig

R = results()

st.subheader("What we learned from one slice of tonsil")
with st.container(horizontal=True):
    for icon, title, text in [
        (":material/map:", "Genes + position", "Spatial transcriptomics measures thousands of genes *and* where they are."),
        (":material/cleaning_services:", "Clean & adjust", "Remove poor spots and make spots comparable."),
        (":material/bubble_chart:", "Genes alone draw anatomy", "Grouping spots by genes only gives real tissue regions."),
        (":material/label:", "Named & confirmed", "Marker genes name each region; proteins confirm it."),
        (":material/hub:", "Organisation", "Germinal centres sit inside follicles; T cells keep away from the surface."),
    ]:
        with st.container(border=True, width=210, height="stretch"):
            st.markdown(f"### {icon}")
            st.markdown(f"**{title}**")
            st.caption(text)

st.space("medium")
c1, c2 = st.columns([1.2, 1])
with c1:
    fig = tissue_fig(height=480)
    add_region_traces(fig)
    fig.update_layout(legend=dict(orientation="h", y=-0.02))
    st.plotly_chart(fig, alt="Final named tissue map")
with c2:
    st.markdown("**:material/warning: Keep in mind**  \n"
                "One spot = several cells, and this is one sample from one person.")
    st.markdown("**:material/arrow_forward: Next in the field**  \n"
                "Single-cell resolution (Visium HD, Xenium, MERFISH) and estimating cell types per spot.")
    st.space("small")
    references("stahl2016", "marx2021", "massoni2024", "wolf2018", "palla2022", "dataset")
