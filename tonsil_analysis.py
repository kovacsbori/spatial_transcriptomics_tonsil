"""Analysis pipeline for the 10x CytAssist FFPE human tonsil dataset (RNA + protein).

What this file does, in order:
    1. Load the RNA counts, the spot positions/H&E image and the antibody (protein) counts.
    2. Quality control: measure every spot, drop the ones with too little signal.
    3. Normalise RNA, pick informative genes, compress (PCA), build a similarity graph, cluster (Leiden)
       and draw a 2-D map (UMAP).
    4. Name the clusters by scoring known tonsil marker-gene sets.
    5. Normalise the protein counts and compare them with RNA from the same spots.
    6. Spatial statistics: which regions touch (neighbourhood enrichment) and which genes form
       patches (Moran's I).

Every intermediate result is kept (raw counts, QC before filtering, graphs, scores) so the Streamlit app
can show each step. Run directly (`python tonsil_analysis.py`) to print the key results.

Data layout reminder (AnnData):
    adata.X        spots × genes matrix (rows = spots, columns = genes)
    adata.obs      one row per spot: QC metrics, cluster labels, scores …
    adata.var      one row per gene: gene metadata, "highly variable" flag …
    adata.obsm     per-spot arrays: "spatial" (x, y in full-resolution pixels), "X_pca", "X_umap"
    adata.obsp     spot × spot graphs: expression neighbours, spatial neighbours
    adata.uns      everything else: images, scale factors, results of statistics
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import squidpy as sq
from scipy import sparse
from scipy.stats import spearmanr

# scanpy/squidpy/anndata print many FutureWarnings about upcoming API changes. They don't affect the
# results, and they would clutter the app's terminal output, so we silence them.
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------------------------------- settings
# Folder with Space Ranger's output: filtered_feature_bc_matrix.h5 + spatial/ (image, positions, scale factors)
DATA_DIR = Path(__file__).parent / "data" / "Human_Tonsil_FFPE"

# Fixed random seed: PCA, UMAP, Leiden and gene scoring all involve randomness. Fixing the seed makes every
# run give exactly the same clusters and pictures (important: the app's text refers to specific clusters).
SEED = 0

# QC thresholds: a spot must have at least this many RNA molecules (UMIs) and distinct genes to be kept.
# They are deliberately mild: in spatial data, low RNA is often real biology (e.g. connective tissue has
# fewer cells per spot), not a technical failure.
QC_MIN_COUNTS, QC_MIN_GENES = 500, 250

# Leiden "resolution": higher = more, smaller clusters. 0.3 gives 6 clusters (one per named tissue part).
# Note: at 0.6 (8 clusters) the compartments separate more cleanly; at 6 the T-cell zone and its stroma
# stay together and plasma cells mix with epithelium.
LEIDEN_RESOLUTION = 0.3

# Known tonsil programs (textbook marker genes) used to give clusters a name. Each list should contain genes
# that are specifically high in that compartment. Sources: tonsil atlas (Massoni-Badosa et al. 2024) and
# standard immunology markers.
SIGNATURES = {
    # resting / naive B cells forming the outer ring of a follicle
    "B follicle (mantle)": ["MS4A1", "CD19", "CD79A", "FCER2", "IGHD", "TCL1A", "IGHM"],
    # GC B-cell program (BCL6, AICDA, RGS13, LMO2, MEF2B), proliferation (MKI67, TOP2A) and
    # follicular dendritic cells that build the GC scaffold (CR2, FDCSP, CXCL13)
    "Germinal center": ["BCL6", "AICDA", "RGS13", "LMO2", "MEF2B", "MKI67", "TOP2A", "CR2", "FDCSP", "CXCL13"],
    # T cells (CD3E, CD3D, CD4, IL7R) and the chemokines that hold them in the T-cell zone (CCL19, CCL21 -> CCR7)
    "T-cell zone": ["CD3E", "CD3D", "CD4", "CCR7", "IL7R", "CCL19", "CCL21", "LTB"],
    # antibody-secreting cells: J chain, ER/secretion machinery (MZB1, XBP1, PRDM1) and antibody genes
    "Plasma cells": ["JCHAIN", "MZB1", "XBP1", "PRDM1", "IGHG1", "IGKC"],
    # stratified squamous epithelium lining the tonsil surface and crypts (keratins, cornification genes)
    "Squamous epithelium": ["KRT5", "KRT13", "KRT14", "KRT4", "SPRR3", "S100A2"],
    # connective tissue (collagens), smooth muscle (ACTA2, MYH11) and blood vessels (PECAM1, VWF)
    "Stroma / vessels": ["COL1A1", "COL3A1", "ACTA2", "MYH11", "PECAM1", "VWF"],
}


# ------------------------------------------------------------------------------------------- loading
def load():
    """Return (adata, prot): RNA with spatial info, and protein (antibody capture) for the same spots."""
    # read_visium reads the gene-expression part of the .h5 file and attaches the spatial folder:
    #   adata.obsm["spatial"]  -> spot centres in full-resolution image pixels
    #   adata.uns["spatial"]   -> the downscaled H&E images and the scale factors to map between them
    adata = sc.read_visium(DATA_DIR)
    # A few gene symbols appear twice in the probe panel; make them unique (GENE, GENE-1, …)
    adata.var_names_make_unique()

    # The same .h5 file also contains the antibody counts ("Antibody Capture"). read_visium drops them,
    # so we read the file again with gex_only=False and keep only the antibody columns.
    full = sc.read_10x_h5(DATA_DIR / "filtered_feature_bc_matrix.h5", gex_only=False)
    prot = full[:, full.var["feature_types"] == "Antibody Capture"].copy()

    # Antibodies are stored under their gene ID (e.g. MS4A1). Keep that for the RNA comparison, but use the
    # better-known protein name (e.g. CD20) as the column name.
    prot.var["gene"] = prot.var["gene_ids"].astype(str)
    prot.var_names = prot.var["secondary_name"].astype(str)
    # Isotype controls are antibodies with no human target: they only measure background binding.
    prot.var["isotype"] = prot.var["isotype_control"].astype(str).str.upper() == "TRUE"
    return adata, prot


# ------------------------------------------------------------------------------------------- pipeline
def run():
    """Run the full pipeline and return everything the app needs in one dictionary."""
    adata, prot = load()

    # The spatial metadata is stored per "library" (one capture area). This dataset has exactly one.
    lib = list(adata.uns["spatial"])[0]
    sp = adata.uns["spatial"][lib]
    # Spot coordinates are in full-resolution pixels; the bundled "hires" image is downscaled.
    # Multiplying coordinates by this factor puts spots on top of the hires image.
    hires_scale = sp["scalefactors"]["tissue_hires_scalef"]

    # === 1. Quality control ==================================================================================
    # Flag mitochondrial genes (names start with "MT-"); a high mitochondrial fraction can indicate damaged
    # tissue. In probe-based FFPE data it is low everywhere, but we still record it.
    adata.var["mt"] = adata.var_names.str.startswith("MT-")
    # Adds per-spot columns to adata.obs: total_counts (UMIs), n_genes_by_counts, pct_counts_mt, …
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Keep a copy of QC for *all* spots (before filtering), with positions, so the app can show
    # which spots a given threshold would remove.
    qc_all = adata.obs[["total_counts", "n_genes_by_counts", "pct_counts_mt"]].copy()
    qc_all[["x", "y"]] = adata.obsm["spatial"] * hires_scale
    n_genes_before = adata.n_vars

    # Remove low-information spots …
    adata = adata[(adata.obs["total_counts"] >= QC_MIN_COUNTS)
                  & (adata.obs["n_genes_by_counts"] >= QC_MIN_GENES)].copy()
    # … and genes detected in fewer than 10 spots (too rare to help find patterns).
    sc.pp.filter_genes(adata, min_cells=10)
    # Keep the protein table in sync: same spots, same order.
    prot = prot[adata.obs_names].copy()

    # === 2. RNA: normalise -> log -> variable genes -> PCA -> graph -> UMAP + Leiden ===========================
    # Keep raw counts in a layer: the HVG method below and the app's "before normalisation" view need them.
    adata.layers["counts"] = adata.X.copy()

    # Normalisation: scale every spot to 10,000 total counts, so spots that captured more RNA overall
    # don't look "higher" for every gene.
    sc.pp.normalize_total(adata, target_sum=1e4)
    # log(1 + x): compresses the range so a few very abundant genes (e.g. immunoglobulins) don't dominate,
    # and makes fold changes comparable at low and high expression. The +1 keeps zeros at zero.
    sc.pp.log1p(adata)

    # Highly variable genes (HVGs): the 2,000 genes that vary more across spots than expected for their
    # mean expression. "seurat_v3" models the mean–variance relationship on raw counts (hence layer="counts").
    # Genes are only flagged in adata.var["highly_variable"]; all genes stay in the object for plotting.
    sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=2000, layer="counts")

    # PCA: compress the 2,000 HVGs into 30 principal components (summary axes of co-varying genes).
    # This keeps the main biological signal and drops much of the noise. Stored in adata.obsm["X_pca"].
    sc.tl.pca(adata, n_comps=30, mask_var="highly_variable", random_state=SEED)

    # Expression neighbour graph: connect each spot to the 15 spots with the most similar PCA profile.
    # This is similarity in *expression*, not physical distance. Stored in adata.obsp["connectivities"].
    sc.pp.neighbors(adata, n_neighbors=15, random_state=SEED)

    # UMAP: 2-D layout of that neighbour graph for plotting (the "similarity map"). Only local closeness
    # is meaningful; axes and distances between far-apart groups are not.
    sc.tl.umap(adata, random_state=SEED)

    # Leiden: find communities of densely connected spots in the neighbour graph = clusters.
    # flavor="igraph" + n_iterations=2 is the fast implementation recommended by current scanpy.
    sc.tl.leiden(adata, resolution=LEIDEN_RESOLUTION, flavor="igraph", n_iterations=2, directed=False,
                 random_state=SEED, key_added="leiden")

    # === 3. Annotation: name clusters with marker-gene signatures =============================================
    # score_genes: per spot, mean expression of the gene set minus the mean of a random reference set of
    # genes with similar expression levels (Tirosh et al. 2016). > 0 = program more active than background.
    # Each score becomes a column in adata.obs named after the signature.
    for name, genes in SIGNATURES.items():
        sc.tl.score_genes(adata, [g for g in genes if g in adata.var_names], score_name=name, random_state=SEED)

    # Average each signature score per cluster, then z-score each signature across clusters so that all
    # signatures are on the same scale (some gene sets give larger raw scores than others).
    cl = adata.obs.groupby("leiden", observed=True)[list(SIGNATURES)].mean()
    cluster_z = (cl - cl.mean()) / cl.std()
    # Each cluster gets the name of its highest-scoring signature; every spot inherits its cluster's name.
    adata.obs["region"] = adata.obs["leiden"].map(cluster_z.idxmax(axis=1)).astype("category")

    # Independent check of the names: the genes that are most specifically higher in each cluster
    # (Wilcoxon rank-sum test, each cluster vs all other spots). use_raw=False -> use the log-normalised X.
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", use_raw=False)
    # Top 6 marker genes per cluster as a table (columns = clusters, rows = rank).
    markers = pd.DataFrame({c: sc.get.rank_genes_groups_df(adata, group=c)["names"].head(6).values
                            for c in adata.obs["leiden"].cat.categories})

    # === 4. Protein: normalise and compare with RNA ===========================================================
    # Dense spots × antibodies matrix of raw counts (small: ~4,200 × 35).
    Xp = prot.X.toarray() if sparse.issparse(prot.X) else np.asarray(prot.X)
    prot_counts = pd.DataFrame(Xp, index=prot.obs_names, columns=prot.var_names)

    # CLR (centred log-ratio) per antibody, the standard for antibody-derived tags (Stoeckius et al. 2017):
    # log(1 + count) minus that antibody's mean over all spots. Each antibody has its own background level;
    # after CLR, 0 = "typical level for this antibody", > 0 = enriched in this spot.
    logx = np.log1p(Xp)
    prot_clr = pd.DataFrame(logx - logx.mean(axis=0, keepdims=True), index=prot.obs_names, columns=prot.var_names)

    # Separate the isotype controls (background only) from antibodies with a real target.
    isotypes = list(prot.var_names[prot.var["isotype"]])
    prot_genes = prot.var["gene"][~prot.var["isotype"]]          # protein name -> gene symbol
    prot_df = prot_clr.drop(columns=isotypes)

    # RNA–protein agreement: for every antibody whose gene is also on the RNA panel, Spearman correlation
    # across all spots between log-normalised RNA and CLR protein. Spearman only assumes a monotonic
    # relationship (higher RNA -> higher protein), not a straight line.
    rows = []
    for p, g in prot_genes.items():
        if g in adata.var_names:
            rna = sc.get.obs_df(adata, keys=[g])[g].values
            rows.append({"protein": p, "gene": g, "rho": spearmanr(rna, prot_df[p].values)[0]})
    rna_prot = pd.DataFrame(rows).sort_values("rho", ascending=False).reset_index(drop=True)

    # === 5. Spatial statistics ===============================================================================
    # Physical neighbour graph: Visium spots sit on a hexagonal grid, so each spot has up to 6 direct
    # neighbours. Stored in adata.obsp["spatial_connectivities"] (separate from the expression graph).
    sq.gr.spatial_neighbors(adata, coord_type="grid", n_neighs=6)

    # Neighbourhood enrichment: count edges between each pair of regions in the spatial graph and compare with
    # 1,000 random shuffles of the region labels. Result is a z-score per pair:
    # > 0 = regions touch more than chance, < 0 = they avoid each other.
    sq.gr.nhood_enrichment(adata, cluster_key="region", seed=SEED, show_progress_bar=False)
    cats = adata.obs["region"].cat.categories
    nhood = pd.DataFrame(adata.uns["region_nhood_enrichment"]["zscore"], index=cats, columns=cats)

    # Moran's I for every HVG: does a spot's value resemble its spatial neighbours' values?
    # ~1 = smooth patches, ~0 = random. n_perms=None uses the fast analytic p-value instead of permutations.
    # Results are sorted from most to least spatially patterned.
    hvg = adata.var_names[adata.var["highly_variable"]].tolist()
    sq.gr.spatial_autocorr(adata, mode="moran", genes=hvg, n_perms=None, n_jobs=1)
    moran = adata.uns["moranI"][["I"]].copy()

    # === 6. Image and coordinates for plotting ===============================================================
    # The hires H&E image is stored as floats in [0, 1]; convert to 8-bit for fast JPEG encoding in Plotly.
    image = (np.clip(sp["images"]["hires"], 0, 1) * 255).astype(np.uint8)
    # Spot centres in hires-image pixels (so they line up with `image`).
    xy = adata.obsm["spatial"] * hires_scale

    return dict(
        adata=adata,                    # filtered, normalised RNA with clusters, scores, graphs, UMAP
        prot=prot_df,                   # CLR-normalised protein (isotype controls removed)
        prot_counts=prot_counts,        # raw antibody counts incl. isotype controls
        isotypes=isotypes,              # names of the isotype-control antibodies
        qc_all=qc_all,                  # QC metrics + positions for ALL spots (before filtering)
        n_genes_before=n_genes_before,  # number of genes before gene filtering
        cluster_z=cluster_z,            # cluster × signature z-scored mean score (used for naming)
        markers=markers,                # top Wilcoxon marker genes per cluster
        rna_prot=rna_prot,              # RNA–protein Spearman correlation per antibody
        nhood=nhood,                    # region × region neighbourhood-enrichment z-scores
        moran=moran,                    # Moran's I per HVG, most patterned first
        image=image,                    # hires H&E image (uint8)
        xy=xy,                          # spot centres in hires-image pixels
        spot_diameter_px=sp["scalefactors"]["spot_diameter_fullres"] * hires_scale,  # spot size on the image
    )


if __name__ == "__main__":
    # Quick text summary of the main results, useful to check the pipeline without starting the app.
    r = run()
    a = r["adata"]
    print(a)
    print("\nSpots per region:\n", a.obs["region"].value_counts().to_string())
    print("\nMarkers:\n", r["markers"].to_string())
    print("\nRNA vs protein:\n", r["rna_prot"].round(2).head(10).to_string())
    print("\nNhood z:\n", r["nhood"].round(1).to_string())
    print("\nQC:\n", r["qc_all"].describe().round(1).to_string())
