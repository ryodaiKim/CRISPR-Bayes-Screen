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

How the pipeline fits together (quick workflow)
----------------------------------------------

- 1) Data: place raw input CSVs in `data/` (e.g. `sgrna_summary.csv`, `gene_summary.csv`).
- 2) Preprocessing: run the notebook cells that build `counts_long` (tidy counts), compute library-size `size_factor`, and save `outputs/tables/tidy_counts.csv`.
- 3) Empirical-Bayes (EB) calibration: run the EB calibration cell to estimate global hyperparameter scales (s_beta, s_u, sigma_logphi) and a LOWESS smoother for gene-wise phi. This writes `outputs/tables/global_params_priors.json` and `outputs/tables/mu_logphi_loess.joblib`.
- 4) Fast stage-1 screening (candidate selection): run the candidate generation cells to compute per-gene metrics and `composite_rank`, producing `outputs/tables/candidate_genes.csv` and `outputs/tables/tidy_counts_candidates.csv`.
- 5) Stratified batch assignment (for joint NUTS): run the stratified assignment cell which divides the top 6,000 candidate genes into 10 deciles by `composite_rank` and constructs 10 batches of 600 genes each (60 genes per decile per batch). This writes `outputs/tables/batch_assignments_stratified_by_composite_rank.csv`.
- 6) Joint-batch NUTS sampling: run the NUTS cell which accepts an optional `BATCH_NO` (1..10). If `BATCH_NO` is provided, it reads the stratified assignment CSV and jointly samples that batch under the HBM with EB-informed global priors. Per-batch outputs:
	- `outputs/tables/posterior_cache_hbm_nuts/beta_draws_batch_{BATCH_NO:04d}.npz` (posterior draws cache),
	- updates `outputs/tables/bayes_gene_summary_hbm_nuts.csv` (per-gene summary),
	- accumulates idata in `outputs/tables/posterior_cache_hbm_nuts/idata_gene_accum.nc`.
- 7) Diagnostics & hit calling: run diagnostics cells (R-hat, ESS) and the hit-calling cells which aggregate per-gene probabilities and write `hits_hbm_nuts.csv`, augmented summaries, and volcano/diagnostic figures under `outputs/`.

Notes & tips
------------
- The stratified batching uses deterministic sampling (fixed random seed) so assignments are reproducible. If you modify `TOTAL_GENES`, `BATCHES` or `PER_DECILE_PER_BATCH`, re-run the assignment cell to regenerate new batches.
- NUTS sampling is computationally intensive; consider running one batch interactively first to validate configuration (chains, draws, tune, target_accept), then scale to multiple batches.
- Outputs are resume-safe: `bayes_gene_summary_hbm_nuts.csv` is merged and `idata_gene_accum.nc` is appended/updated so re-running batches won't re-sample genes already present unless `OVERWRITE_EXISTING=True` in the NUTS cell.

If you'd like, I can also add a small `scripts/` runner to invoke the notebook cells from the command line (nbclient/nbconvert) to run batches headlessly.