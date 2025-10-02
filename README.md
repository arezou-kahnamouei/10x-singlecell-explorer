# 10x Single-cell Explorer (Seurat + Shiny)

QC → Clustering → UMAP, VPS36 stats, per-cluster DE, gene-set DE, advanced violins.

## Quickstart
```bash
# Create env (optional)
conda env create -f environment.yml
conda activate seurat_app

# Run locally
R -q -e "shiny::runApp('.', launch.browser=TRUE)"
