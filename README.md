# CRISPR-Bayes-Screen

This repository implements a Hierarchical Bayesian Model (HBM) for pooled CRISPR screens and supporting benchmark pipelines. The primary goal is joint modeling of all genes and sgRNAs under a Negative-Binomial observation model with global shrinkage on per-gene effects. The project contains both variational (VI) and full MCMC (NUTS) workflows, utilities for fast stage‑1 screening, and comparison pipelines including an empirical-Bayes ShrinkCRISPR-style implementation and a MAGeCK benchmark.

Key ideas and outputs
- HBM: a hierarchical NB-GLM with gene intercepts, per-gene treatment effects \u03b2_g, sgRNA random effects, and sample offsets (library-size factors). Inference is done with ADVI (fast, cached draws) and optional NUTS refinement for a selected gene set. Outputs include per-gene summaries, posterior caches (NPZ), and hit lists based on lfsr and probability-of-effect thresholds.
- ShrinkCRISPR (EB benchmark) *yet to be developed: a standalone R pipeline implementing ShrinkBayes/INLA-style Bayes factors → lfdr decisions. Because INLA/fmesher builds must match exactly, an in-notebook rpy2→INLA flow may fail in some environments; a fallback R script is provided (see `scripts/shrinkcrispr_inla.R`).
- MAGeCK benchmark: MAGeCK gene-summary conversion and hit tables are included for comparison.

Repository layout
- `data/` — input CSVs: `sgrna_summary.csv`, `gene_summary.csv`, reference lists.
- `notebooks/CRISPR_Bayes/` — analysis notebooks (HBM pipeline, GLM analysis, ShrinkCRISPR attempts).
- `scripts/` — helper scripts (including `scripts/shrinkcrispr_inla.R`, an R fallback that runs INLA externally and writes NPZ/CSV outputs). *yet to be developed
- `outputs/figures/`, `outputs/tables/` — generated results and caches (posterior caches are under `outputs/tables/posterior_cache_*`).
- `requirements.txt` — Python dependencies for notebook runs.

Outputs and hit-calling
- Per-gene VI summary: `outputs/tables/bayes_gene_summary_hbm_vi.csv` (and augmented versions ending in `_with_probs.csv`).
- Per-gene NUTS summary: `outputs/tables/bayes_gene_summary_hbm_nuts.csv` and `_with_probs.csv` plus `posterior_cache_hbm_nuts/` NPZ files.
- Hit lists are produced as `hits_*.csv` with configurable thresholds (lfsr and P(|beta| >= delta)). The default notebook thresholds are documented in the notebooks; they are intentionally configurable.