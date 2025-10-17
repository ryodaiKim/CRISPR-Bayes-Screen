# CRISPR-Bayes-Screen

This repo centers on a Negative Binomial (NB) multi-gene Bayesian model for CRISPR screens. A ShrinkCRISPR-aligned empirical-Bayes pipeline is included as a benchmark for accuracy and calibration.

What’s here:
- HBM Bayesian model (main focus): joint inference across genes/sgRNAs with VI/NUTS, caching, robust hit-calling, QC and volcano plots.
- ShrinkCRISPR benchmark (supporting): empirical-Bayes approximation aligned with the paper (zero-centered prior, lfdr-based calls) used to compare against the NB model.
 - MAGeCK benchmark (supporting): a commonly used, alternative gene-level calling method used in our comparisons and benchmarking (see `outputs/tables/bench_mageck_hits.csv` and related bench outputs).

## Layout
- `data/`: input CSVs (sgRNA/gene summaries, reference sets)
- `notebooks/CRISPR_Bayes/`: analysis notebooks (sampling, GLM analysis, ShrinkCRISPR)
- `outputs/figures/` and `outputs/tables/`: generated plots and tables
- `src/`: future library code (if consolidated)

## Quick start
Run the notebooks in `notebooks/CRISPR_Bayes/`:

1) `01_bayesian_screen.ipynb` — Fit the HBM model (VI/NUTS), save posterior caches and summaries.
2) `02_nb_glm_analysis.ipynb` — Analyze/visualize NB results: QC, global histograms, volcano plots, and overlaps.
3) `03_ShrinkCRISPR.ipynb` — Run the ShrinkCRISPR-aligned EB pipeline to benchmark the NB model.

Artifacts are saved under `outputs/` automatically.

## Environment
Install Python dependencies:

- See `requirements.txt`

If using notebooks, ensure the kernel matches the project environment.

## Benchmarks and thresholds
- NB model hit-calling: based on posterior-based lfsr and probability of effect (configurable in notebooks), with non-hits shown in gray on volcano plots.
- ShrinkCRISPR benchmark: EB with zero-centered prior and lfdr-based decisions (BH q-values as a practical proxy). Default lfdr threshold: 0.10 (paper used 0.05 in sims, 0.10 in examples). Optional probability-of-effect filter can be enabled for conservativeness.
 - MAGeCK: included as an external benchmark for gene calling; many of the benchmarking figures and overlap tables compare NB, ShrinkCRISPR and MAGeCK results.

## Notes
- Large posterior caches and checkpoints are ignored via `.gitignore`.
- For rscreenorm-style lethality mapping in ShrinkCRISPR, provide negative/positive control gene lists and enable the option in the notebook. If T0 counts are available, full fold-change + rscreenorm can be enabled.