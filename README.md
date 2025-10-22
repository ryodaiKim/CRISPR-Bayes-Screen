# CRISPR-Bayes-Screen

This repository implements a Hierarchical Bayesian Model (HBM) for pooled CRISPR screens and supporting benchmark pipelines. The primary goal is joint modeling of all genes and sgRNAs under a Negative-Binomial observation model with global shrinkage on per-gene effects. 

Key ideas and outputs
- HBM: a hierarchical NB-GLM with gene intercepts, per-gene treatment effects $β_g$, sgRNA random effects, and sample offsets (library-size factors). Inference is done with NUTS for a selected gene set. Outputs include per-gene summaries, posterior caches (NPZ), and hit lists based on lfsr and probability-of-effect thresholds.
- ShrinkCRISPR (EB benchmark) *yet to be developed: a standalone R pipeline implementing ShrinkBayes/INLA-style Bayes factors → lfdr decisions. Because INLA/fmesher builds must match exactly, an in-notebook rpy2→INLA flow may fail in some environments; a fallback R script is provided (see `scripts/shrinkcrispr_inla.R`).
- MAGeCK benchmark: MAGeCK gene-summary conversion and hit tables are included for comparison.

- NEL (Laplace approximation) updates: the per-gene NB-GLM MAP solver now supports a configurable prior on the treatment effect β. Choose between
  - Gaussian prior (ridge): `β ~ Normal(0, τ_beta^2)` (previous default), and
  - Laplace prior (L1 via LQA): `β ~ Laplace(0, b)` with `Var=2 b^2` and `b = τ_beta / √2`.
  The Laplace option yields heavier tails and moderates shrinkage of large effects. It is enabled in the notebook call to `fit_all_genes(..., prior_beta='laplace')`. Set `prior_beta='gaussian'` to revert.

Repository layout
- `data/` — input CSVs: `sgrna_summary.csv`, `gene_summary.csv`, reference lists.
- `notebooks/CRISPR_Bayes/` — analysis notebooks (HBM pipeline, GLM analysis, ShrinkCRISPR attempts).
- `scripts/` *yet to be developed — helper scripts (including `scripts/shrinkcrispr_inla.R`, an R fallback that runs INLA externally and writes NPZ/CSV outputs).
- `outputs/figures/`, `outputs/tables/` — generated results and caches (posterior caches are under `outputs/tables/posterior_cache_*`).
- `requirements.txt` — Python dependencies for notebook runs.

Outputs and hit-calling
- Per-gene NUTS summary: `outputs/tables/bayes_gene_summary_hbm_nuts.csv` and `_with_probs.csv` plus `posterior_cache_hbm_nuts/` NPZ files.
- Hit lists are produced as `hits_*.csv` with configurable thresholds (lfsr and P(|beta| >= delta)). The default notebook thresholds are documented in the notebooks; they are intentionally configurable.

NEL shrinkage control (β prior)
- Tune `τ_beta` to control shrinkage strength for β in the NEL solver. Larger `τ_beta` → less shrinkage; smaller `τ_beta` → more.
- With `prior_beta='laplace'`, the solver uses a stable local quadratic approximation (LQA) to implement an L1 penalty on β. sgRNA effects `u` remain Gaussian (ridge) with scale `τ_u`.
- You can override the empirical estimate: e.g. `priors['tau_beta'] = 0.6` or multiply it (e.g. `* 1.5`) before calling `fit_all_genes`.

How the pipeline fits together (quick workflow)
----------------------------------------------

- 1) Data: place raw input CSVs in `data/` (e.g. `sgrna_summary.csv`, `gene_summary.csv`).
- 2) Preprocessing: run the notebook cells that build `counts_long` (tidy counts), compute library-size `size_factor`, and save `outputs/tables/tidy_counts.csv`.
- 3) Empirical-Bayes (EB) calibration: run the EB calibration cell to estimate global hyperparameter scales (s_beta, s_u, sigma_logphi) and a LOWESS smoother for gene-wise phi. This writes `outputs/tables/global_params_priors.json` and `outputs/tables/mu_logphi_loess.joblib`.
- 4) Stratified batch assignment (for joint NUTS): run the stratified assignment cell which divides the top 6,000 candidate genes into 10 deciles by `composite_rank` and constructs 10 batches of 600 genes each (60 genes per decile per batch). This writes `outputs/tables/batch_assignments_stratified_by_composite_rank.csv`.
- 5) Joint-batch NUTS sampling: run the NUTS cell which accepts an optional `BATCH_NO` (1..10). If `BATCH_NO` is provided, it reads the stratified assignment CSV and jointly samples that batch under the HBM with EB-informed global priors. Per-batch outputs:
	- `outputs/tables/posterior_cache_hbm_nuts/beta_draws_batch_{BATCH_NO:04d}.npz` (posterior draws cache),
	- updates `outputs/tables/bayes_gene_summary_hbm_nuts.csv` (per-gene summary),
	- accumulates idata in `outputs/tables/posterior_cache_hbm_nuts/idata_gene_accum.nc`.
- 6) Diagnostics & hit calling: run diagnostics cells (R-hat, ESS) and the hit-calling cells which aggregate per-gene probabilities and write `hits_hbm_nuts.csv`, augmented summaries, and volcano/diagnostic figures under `outputs/`.

Notes & tips
------------
- The stratified batching uses deterministic sampling (fixed random seed) so assignments are reproducible. If you modify `TOTAL_GENES`, `BATCHES` or `PER_DECILE_PER_BATCH`, re-run the assignment cell to regenerate new batches.
- NUTS sampling is computationally intensive; consider running one batch interactively first to validate configuration (chains, draws, tune, target_accept), then scale to multiple batches.
- Outputs are resume-safe: `bayes_gene_summary_hbm_nuts.csv` is merged and `idata_gene_accum.nc` is appended/updated so re-running batches won't re-sample genes already present unless `OVERWRITE_EXISTING=True` in the NUTS cell.
