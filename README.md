# 10x Single-cell Explorer (Seurat + Shiny)

QC → Clustering → UMAP, VPS36 stats, per-cluster DE, gene-set DE, advanced violins.

## Quickstart
```bash
# Create env (optional)
conda env create -f environment.yml
conda activate seurat_app

# Run locally
R -q -e "shiny::runApp('.', launch.browser=TRUE)"
# 10x Single-cell Explorer (Seurat + Shiny)

[![GitHub release (latest by date)](https://img.shields.io/github/v/release/USERNAME/10x-singlecell-explorer)](https://github.com/USERNAME/10x-singlecell-explorer/releases)
[![CI](https://github.com/USERNAME/10x-singlecell-explorer/actions/workflows/shiny-smoke.yml/badge.svg)](https://github.com/USERNAME/10x-singlecell-explorer/actions/workflows/shiny-smoke.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

QC → Clustering → UMAP, VPS36 stats, per-cluster DE, gene-set DE, advanced violins.
