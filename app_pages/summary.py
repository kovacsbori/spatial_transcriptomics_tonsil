import streamlit as st

from tonsil_ui import add_region_traces, references, results, tissue_fig

R = results()

c1, c2, c3 = st.columns([1, 1.1, 1])
with c1:
    st.subheader("What we learned")
    for icon, title, text in [
        (":material/map:", "Genes + position", "Thousands of genes *and* where they are."),
        (":material/cleaning_services:", "Clean & adjust", "Remove poor spots, make spots comparable."),
        (":material/bubble_chart:", "Genes alone draw anatomy", "Grouping by genes only gives real tissue regions."),
        (":material/label:", "Named & confirmed", "Marker genes name each region; proteins confirm it."),
        (":material/hub:", "Organisation", "Germinal centres inside follicles; T cells away from the surface."),
    ]:
        st.markdown(f"{icon} **{title}**  \n:gray[{text}]")
with c2:
    fig = tissue_fig(height=470)
    add_region_traces(fig)
    fig.update_layout(legend=dict(orientation="h", y=-0.02))
    st.plotly_chart(fig, alt="Final named tissue map")
with c3:
    st.markdown("**:material/warning: Keep in mind**  \n"
                ":gray[One spot = several cells, and this is one sample from one person.]")
    st.markdown("**:material/arrow_forward: Next in the field**  \n"
                ":gray[Single-cell resolution (Visium HD, Xenium, MERFISH) and estimating cell types per spot.]")
    references("stahl2016", "marx2021", "massoni2024", "wolf2018", "palla2022", "dataset")
