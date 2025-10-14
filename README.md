# CRISPR-Bayes-Screen

Bayesian and empirical-Bayes pipelines for CRISPR screen analysis:
- Negative Binomial multi-gene model with NUTS/VI and robust hit-calling
- ShrinkCRISPR-style empirical-Bayes approximation with lfdr-based decisions

## Layout
- `data/`: input CSVs (sgRNA/gene summaries, reference sets)
- `notebooks/CRISPR_Bayes/`: analysis notebooks (sampling, GLM analysis, ShrinkCRISPR)
- `outputs/figures/` and `outputs/tables/`: generated plots and tables
- `src/`: future library code (if consolidated)

## Quick start
Open and run the notebooks in `notebooks/CRISPR_Bayes/` in order:
1) `01_bayesian_screen.ipynb` (model fit/caching)
2) `02_nb_glm_analysis.ipynb` (QC, volcano, overlaps)
3) `03_ShrinkCRISPR.ipynb` (EB approximation aligned with ShrinkCRISPR paper)

Artifacts are saved under `outputs/` automatically.

## Environment
Install Python dependencies:

- See `requirements.txt`

If using notebooks, ensure the kernel matches the project environment.

## Notes
- Large posterior caches and checkpoints are ignored via `.gitignore`.
- For rscreenorm-style lethality mapping, provide negative/positive control gene lists and enable the option in the ShrinkCRISPR notebook.