# Spatial transcriptomics walkthrough: human tonsil

An interactive Streamlit presentation that introduces spatial transcriptomics using a public 10x Genomics
Visium CytAssist dataset of human tonsil (RNA + 31 proteins measured on the same tissue section).

It walks through the analysis in nine slides: what spatial transcriptomics is, what a spot contains,
quality control and normalisation, clustering, naming clusters with marker genes, linking expression
space to tissue space, RNA vs protein, and spatial neighbourhoods.

## Setup

Requires Python 3.12. [uv](https://docs.astral.sh/uv/) downloads it automatically if it isn't installed.

```bash
# 1. Install uv (skip if you already have it)
curl -LsSf https://astral.sh/uv/install.sh | sh      # macOS / Linux
#   or: brew install uv                               # macOS with Homebrew
#   or (Windows PowerShell): powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"

# 2. Create an environment and install dependencies
uv venv --python 3.12 .venv
uv pip install --python .venv -r requirements.txt
#   (without uv: python3.12 -m venv .venv && .venv/bin/pip install -r requirements.txt)

# 3. Download the dataset (~70 MB) into data/Human_Tonsil_FFPE/
uv run download_data.py

# 4. Start the app
uv run streamlit run tonsil_app.py
```

The first page load runs the full analysis (~25 seconds); after that every slide is instant.

## Project layout

| Path | Contents |
|---|---|
| `tonsil_app.py` | Entry point: slide order, sidebar, Back/Next navigation |
| `app_pages/` | One file per slide |
| `tonsil_analysis.py` | The analysis pipeline (scanpy + squidpy); run it directly to print key results |
| `tonsil_ui.py` | Shared plotting and layout helpers |
| `tonsil_refs.py` | References shown on the final slide |
| `download_data.py` | Fetches the dataset from 10x Genomics |
| `.streamlit/config.toml` | Theme (colours, fonts) |

## Data

10x Genomics, *CytAssist FFPE Protein Expression – Human Tonsil* (Space Ranger 2.1.0), public dataset:
<https://www.10xgenomics.com/datasets>. The data is not stored in this repository.

## Tools

[scanpy](https://doi.org/10.1186/s13059-017-1382-0) · [squidpy](https://doi.org/10.1038/s41592-021-01358-2) ·
[Streamlit](https://streamlit.io) · [Plotly](https://plotly.com/python/)
