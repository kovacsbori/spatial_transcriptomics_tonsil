import plotly.graph_objects as go
import streamlit as st

from tonsil_ui import REGION_COLORS, how_it_works, in_crop, question, region_crop, results, tissue_fig

R = results()

question("Which genes are switched on, and *where* in the tissue?")

st.markdown("Think of three ways to study a fruit tart:")
with st.container(horizontal=True):
    for icon, name, text in [
        (":material/blender:", "Bulk RNA-seq = a smoothie",
         "Blend the whole tissue. You know *what* is in it, but everything is mixed together."),
        (":material/nutrition:", "Single-cell RNA-seq = a fruit salad",
         "You see every piece, but not *where* it was on the tart."),
        (":material/bakery_dining:", "Spatial transcriptomics = the tart itself",
         "You see what is there **and** where it sits."),
    ]:
        with st.container(border=True, width=300, height="stretch"):
            st.markdown(f"### {icon}")
            st.markdown(f"**{name}**")
            st.caption(text)

st.space("medium")
st.subheader("Our sample: a slice of human tonsil")
c1, c2 = st.columns([1.3, 1])
with c1:
    show_spots = st.toggle("Show the measurement spots", value=True)
    fig = tissue_fig(height=540, force_image=True)
    if show_spots:
        fig.add_trace(go.Scattergl(x=R["xy"][:, 0], y=R["xy"][:, 1], mode="markers", hoverinfo="skip",
                                   marker=dict(size=4, color="rgba(30,30,30,0.45)"), showlegend=False))
    st.plotly_chart(fig, alt="Stained image of the tonsil slice with measurement spots")
with c2:
    zoom = st.segmented_control("Zoom in on", ["Germinal center", "T-cell zone", "Squamous epithelium"],
                                default="Germinal center", required=True)
    crop = region_crop(zoom, half_width=110)
    zfig = tissue_fig(height=420, crop=crop, force_image=True)
    m = in_crop(crop)
    px_per_unit = 400 / (crop[1] - crop[0])
    zfig.add_trace(go.Scatter(x=R["xy"][m, 0], y=R["xy"][m, 1], mode="markers", hoverinfo="skip", showlegend=False,
                              marker=dict(size=R["spot_diameter_px"] * px_per_unit, symbol="circle-open",
                                          line=dict(width=2, color=REGION_COLORS[zoom]))))
    st.plotly_chart(zfig, alt=f"Zoomed view of the {zoom.lower()} with spot outlines")
    st.caption(f"Each circle is one spot, about half the width of a hair, covering **several cells**. "
               f"{len(R['qc_all']):,} spots × ~18,000 genes + {R['prot'].shape[1]} proteins.")

how_it_works([
    "A thin slice of tonsil (5 µm, preserved in paraffin) is put on a glass slide, stained pink/purple (H&E) "
    "and photographed.",
    "Small DNA probes stick to the RNA of ~18,000 genes; DNA-tagged antibodies stick to 31 proteins.",
    "The slice is pressed onto a slide covered with ~5,000 spots. Each spot has its own DNA 'postcode' "
    "(spatial barcode) that is attached to everything it catches.",
    "Sequencing reads out *which gene* and *which postcode*. Counting them gives a table of genes per spot, "
    "and the postcodes put each spot back onto the photo.",
    "Technology: 10x Genomics Visium CytAssist. Public dataset: 'CytAssist FFPE Protein Expression – Human Tonsil'.",
])
