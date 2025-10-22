# CRISPR-Bayes-Screen

This repository implements a Hierarchical Bayesian Model (HBM) for pooled CRISPR screens and supporting benchmark pipelines. The primary goal is joint modeling of all genes and sgRNAs under a Negative-Binomial observation model with global shrinkage on per-gene effects. 

Key ideas and outputs
- HBM: a hierarchical NB-GLM with gene intercepts, per-gene treatment effects $β_g$, sgRNA random effects, and sample offsets (library-size factors). Inference is done with NUTS for a selected gene set. Outputs include per-gene summaries, posterior caches (NPZ), and hit lists based on lfsr and probability-of-effect thresholds.

- NEL (Laplace approximation) updates: the per-gene NB-GLM MAP solver now supports a configurable prior on the treatment effect β. Choose between
  - Gaussian prior (ridge): `β ~ Normal(0, τ_beta^2)` (previous default), and
  - Laplace prior (L1 via LQA): `β ~ Laplace(0, b)` with `Var=2 b^2` and `b = τ_beta / √2`.
  The Laplace option yields heavier tails and moderates shrinkage of large effects. It is enabled in the notebook call to `fit_all_genes(..., prior_beta='laplace')`. Set `prior_beta='gaussian'` to revert.

- ShrinkCRISPR (EB benchmark) *yet to be developed: a standalone R pipeline implementing ShrinkBayes/INLA-style Bayes factors → lfdr decisions. Because INLA/fmesher builds must match exactly, an in-notebook rpy2→INLA flow may fail in some environments; a fallback R script is provided (see `scripts/shrinkcrispr_inla.R`).

- MAGeCK benchmark: MAGeCK gene-summary conversion and hit tables are included for comparison.

Repository layout
- `data/` — input CSVs: `sgrna_summary.csv`, `gene_summary.csv`, reference lists.
- `notebooks/CRISPR_Bayes/` — analysis notebooks (HBM and NEL pipeline, ShrinkCRISPR attempts).
- `outputs/figures/`, `outputs/tables/` — generated results and caches (posterior caches are under `outputs/tables/posterior_cache_*`).
- `requirements.txt` — Python dependencies for notebook runs.