import plotly.graph_objects as go
import streamlit as st

from tonsil_ui import REGION_COLORS, how_it_works, in_crop, question, region_crop, results, tissue_fig

R = results()

question("Which genes are switched on, and *where* in the tissue?")

c1, c2, c3 = st.columns([1, 1.15, 1.05], gap="medium")
with c1:
    # Three ways to study a fruit tart, as compact rows
    for icon, name, text in [
        (":material/blender:", "Bulk RNA-seq = a smoothie", "You know *what* is in it, but everything is mixed."),
        (":material/nutrition:", "Single-cell RNA-seq = a fruit salad", "Every piece, but not *where* it was."),
        (":material/bakery_dining:", "Spatial transcriptomics = the tart", "What is there **and** where it sits."),
    ]:
        st.markdown(f"{icon} **{name}**  \n:gray[{text}]")
    how_it_works([
        "5 µm FFPE section, H&E-stained and imaged (10x Visium CytAssist).",
        "Probe pairs hybridise to ~18,000 transcripts; 31 DNA-barcoded antibodies bind proteins.",
        "Probes are captured on ~5,000 spots, each with a unique spatial barcode. Spot **diameter** 55 µm; **distance** between spot centres 100 µm, so the ~45 µm gaps between spots are not measured.",
        "Sequencing reads gene + barcode → a spot × gene count matrix mapped onto the image.",
    ])
with c2:
    show_spots = st.toggle("Show the measurement spots", value=True)
    fig = tissue_fig(height=400, force_image=True, full_image=True)
    if show_spots:
        fig.add_trace(go.Scattergl(x=R["xy"][:, 0], y=R["xy"][:, 1], mode="markers", hoverinfo="skip",
                                   marker=dict(size=4, color="rgba(30,30,30,0.45)"), showlegend=False))
    st.plotly_chart(fig, alt="Stained image of the tonsil slice with measurement spots")
with c3:
    zoom = st.segmented_control("Zoom in on", ["Germinal center", "T-cell zone", "Squamous epithelium"],
                                default="Germinal center", required=True, label_visibility="collapsed",
                                format_func=lambda r: {"Germinal center": "Germ. center", "Squamous epithelium": "Epithelium"}.get(r, r))
    crop = region_crop(zoom, half_width=110)
    zfig = tissue_fig(height=330, crop=crop, force_image=True)
    m = in_crop(crop)
    px_per_unit = (330 - 4) / (crop[1] - crop[0])   # square crop: plot height (minus top margin) in screen px
    zfig.add_trace(go.Scatter(x=R["xy"][m, 0], y=R["xy"][m, 1], mode="markers", hoverinfo="skip", showlegend=False,
                              marker=dict(size=R["spot_diameter_px"] * px_per_unit, symbol="circle-open",
                                          line=dict(width=2, color=REGION_COLORS[zoom]))))
    st.plotly_chart(zfig, alt=f"Zoomed view of the {zoom.lower()} with spot outlines")
    st.caption(f"Each circle is one spot: 55 µm **diameter**, 100 µm centre-to-centre **distance**, covering "
               f"**several cells**. "
               f"{len(R['qc_all']):,} spots × ~18,000 genes + {R['prot'].shape[1]} proteins.")
