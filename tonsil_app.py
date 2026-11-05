"""Guided walkthrough: an introduction to spatial transcriptomics using a human tonsil section.

Run:  .venv/bin/streamlit run tonsil_app.py

This is the entry point. It defines the slide order, the sidebar, the slide header and the Back/Next
buttons. The content of each slide lives in app_pages/<name>.py; shared helpers are in tonsil_ui.py,
the analysis itself in tonsil_analysis.py and the look (colours, fonts) in .streamlit/config.toml.
"""

import streamlit as st

from tonsil_ui import results

# Must be the first Streamlit call: browser-tab title/icon, and a wide layout so maps can sit side by side.
st.set_page_config(
    page_title="Spatial transcriptomics walkthrough",
    page_icon=":material/biotech:",
    layout="wide",
)

# Slide order, grouped into sidebar sections. Each st.Page points to one file in app_pages/.
# To add, remove or reorder slides, edit this dictionary only.
SECTIONS = {
    "The idea": [
        st.Page(
            "app_pages/intro.py",
            title="What is spatial transcriptomics?",
            icon=":material/biotech:",
            default=True,
        ),
        st.Page(
            "app_pages/data.py", title="One spot up close", icon=":material/grid_on:"
        ),
    ],
    "The analysis in 3 steps": [
        st.Page(
            "app_pages/qc.py",
            title="Step 1 · Cleaning the data",
            icon=":material/cleaning_services:",
        ),
        st.Page(
            "app_pages/clustering.py",
            title="Step 2 · Grouping similar spots",
            icon=":material/bubble_chart:",
        ),
        st.Page(
            "app_pages/annotation.py",
            title="Step 3 · Naming the groups",
            icon=":material/label:",
        ),
    ],
    "What we learn": [
        st.Page(
            "app_pages/linked.py",
            title="Similar spots sit together",
            icon=":material/sync_alt:",
        ),
        st.Page(
            "app_pages/protein.py",
            title="RNA vs protein",
            icon=":material/compare_arrows:",
        ),
        st.Page(
            "app_pages/neighbors.py",
            title="Who sits next to whom",
            icon=":material/hub:",
        ),
        st.Page("app_pages/summary.py", title="Summary", icon=":material/flag:"),
    ],
}
# Flat list in presentation order, used for the progress bar and Back/Next.
ALL_PAGES = [p for pages in SECTIONS.values() for p in pages]

# Builds the sidebar menu and returns the page the user is currently on.

page = st.navigation(SECTIONS, position="sidebar", expanded=True)

# Display settings shared by all slides. The `key`s store the values in st.session_state, where
# tonsil_ui.display_opts() reads them when drawing the tissue maps.
with st.sidebar:
    with st.expander("Display settings", icon=":material/tune:"):
        st.slider("Spot size", 2, 14, 6, key="spot_size")
        st.slider("Spot opacity", 0.1, 1.0, 0.85, key="opacity")
        st.toggle("Show H&E image", value=True, key="show_he")

# Position of the current slide in the presentation (0-based).
idx = next(i for i, p in enumerate(ALL_PAGES) if p.url_path == page.url_path)

# Slide header: title + thin progress bar
st.progress((idx + 1) / len(ALL_PAGES))
st.title(page.title)

# Run the analysis (first visit only, ~25 s; cached afterwards), then draw the current slide's content.
results()
page.run()

# Back / Next buttons at the bottom of every slide. "distribute" pushes them to opposite sides.
# Hovering a button shows the title of the slide it leads to.
st.space("medium")
with st.container(horizontal=True, horizontal_alignment="distribute"):
    if idx > 0:
        if st.button(
            "Back", icon=":material/arrow_back:", help=ALL_PAGES[idx - 1].title
        ):
            st.switch_page(ALL_PAGES[idx - 1])
    else:
        st.empty()  # placeholder on the first slide so "Next" still sits on the right
    if idx < len(ALL_PAGES) - 1:
        if st.button(
            "Next",
            icon=":material/arrow_forward:",
            icon_position="right",
            help=ALL_PAGES[idx + 1].title,
        ):
            st.switch_page(ALL_PAGES[idx + 1])
