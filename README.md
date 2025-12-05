# PDC Sampler Condensation Theory & Bioinformatics (ASE vs PDC)

**Version**: 1.0  
**Date**: 2025-11-04  
**First Author**: Gao Han (220224443@seu.edu.cn)  
**Corresponding Author**: Liu Quanjun (lqj@seu.edu.cn)  
**Affiliations**:  
- School of Biological Science and Medical Engineering, Southeast University  
- National Key Laboratory of Digital Medical Engineering

---

## Summary

This repository provides **open, citable** code and selected data associated with our manuscript (under review; DOI/URL will be added upon acceptance).  
It contains two code families:

1) **Condensation theory & Monte Carlo** for the **PDC-sampler** (condensate generation, hygroscopic growth).  
2) **Bioinformatics visualization** for **ASE vs PDC** comparisons.  
   Most upstream computations were performed with **QIIME 2** and **NovoMagic** (Novogene cloud). The R scripts here are primarily for **reproducible visualization and statistics** from exported tables.

Each script begins with a standardized header documenting: **File**, **Title**, **Purpose**, **How to run**, **Inputs/Outputs (repo-relative paths)**, and **the corresponding figure ID(s)** in the manuscript.

**Example (excerpt) – `scripts/bioinfo/WilcoxonGnusD.R`:**

=============================================================================
File: WilcoxonGnusD.R
Title: Genus-level Wilcoxon test (ASE vs PDC) with effect sizes
Purpose
Compute two-sample Wilcoxon rank-sum tests between wet-wall cyclone (ASE)
and PDC-sampler groups at the Genus level, and export a table with:
Taxa, avg(ASE), sd(ASE), avg(PDC), sd(PDC), p.value, q.values,
interval lower, interval upper, log2FoldChange, Combine ES.
The result is intended for manual filtering/sorting in Excel, then
visualization by the companion script WilcoxonGnusF.R (see manuscript Fig. 7a).
Input
./data/bioinfo/featureTable.sample.g.relative.txt
(tab-delimited; columns: Taxonomy | <samples...> [| optional Tax_detail])
How to run
R >= 4.0 with packages: tidyverse, readr, effsize
install.packages(c("tidyverse","readr","effsize"))
From repo root: Rscript scripts/bioinfo/WilcoxonGnusD.R
Output
./results/stats/Wilcoxon_Genus_ASE_vs_PDC_results.csv
Notes
- log2FoldChange = log2( mean(ASE)+eps / mean(PDC)+eps )
- Combine ES = Cliff's delta (ASE vs PDC; >0 means ASE > PDC)
- interval lower/upper = 95% bootstrap CI for (PDC_mean - ASE_mean)
Author : Gao Han
Affil. : School of Biological Science and Medical Engineering, Southeast University
Version: 1.0
Date : 2025-11-04
=============================================================================



> All scripts use **repo-relative paths** like `./data/...`, `./results/...`, `./figures/...` to keep the layout portable and journal-friendly.

---

## Software requirements

- **R ≥ 4.0**  
- Common CRAN packages used across scripts (install on demand in your R session):
  `tidyverse`, `readxl`, `readr`, `ggpubr`, `ggsci`, `svglite`, `effsize`, `lme4`  
- Bioconductor (for ANCOM-BC2 workflow): `phyloseq`, `ANCOMBC`

> Each script’s header lists the minimal packages it needs. No global installer is provided in this repository.

---

## Figure mapping (script → inputs → outputs)

| Manuscript figure(s) | Script | Key input(s) (repo-relative) | Main output(s) (repo-relative) |
|---|---|---|---|
| Fig. 2a | `scripts/theory/CondensationVolume_MentroKaro.R` | — (parameters in script) | `./results/theory/PDC_MonteCarlo_cumulative_time_series.csv`, `./results/theory/PDC_MonteCarlo_summary_by_time.csv` |
| Fig. 2d | `scripts/theory/pdc_gf_curve.R` | — (parameters in script) | `./data/theory/PDC_GF_vs_d0_Sweep.xlsx`, `./figures/fig2/GF_vs_d0_by_S.(svg|pdf|tiff)` |
| Fig. 2e | `scripts/theory/pdc_diameter_time_curves.R` | — (parameters in script) | `./data/theory/PDC_diameter_time_traj_S1p60.xlsx`, `./figures/fig2/diameter_time_S1p60.(svg|pdf|tiff)` |
| Fig. 4a–d | `scripts/bioinfo/AlphaDiversity.R` | `./data/bioinfo/Alpha.xlsx` | `./figures/fig4/Alpha_Boxplots_4in1.(pdf|svg|tiff)` |
| Fig. 4e–g | `scripts/bioinfo/BetaDiversity.R` | `./data/bioinfo/bray_curtis_PCoA.xlsx`, `weighted_unifrac_PCoA.xlsx`, `unweighted_unifrac_PCoA.xlsx` | `./figures/fig4/PCoA_3panels_vertical_BIG.(pdf|svg|tiff)` |
| Fig. 5d | `scripts/bioinfo/Shared_Taxson.R` | `./data/bioinfo/featureTable.group.g.relative.txt` (+ optional shared lists) | `./figures/fig5/Pie_Shared_Genus_byCombinedAbundance.(pdf|svg|tiff)` |
| Fig. 7b, Fig. S3 | `scripts/bioinfo/ANCOMBC2.R` | `./data/bioinfo/featureTable.sample.total.absolute.txt` | `./results/stats/ANCOMBC2_PDC_vs_ASE_results_Genus.csv`, `./figures/fig7/ANCOMBC2_LFC_heatmap_PDC_vs_ASE.(pdf|svg|tiff)` |
| Fig. 7a (table) | `scripts/bioinfo/WilcoxonGnusD.R` | `./data/bioinfo/featureTable.sample.g.relative.txt` | `./results/stats/Wilcoxon_Genus_ASE_vs_PDC_results.csv` |
| Fig. 7a (plots) | `scripts/bioinfo/WilcoxonGnusF.R` | curated `./results/stats/Wilcoxon_Genus_ASE_vs_PDC_results.(xlsx|csv)` | `./figures/fig7/Wilcox_Figure_Panel*.(* )` |
| Fig. 7c | `scripts/bioinfo/WilcoxonPhylumF.R` | curated `./results/stats/Wilcoxon_Phylum_ASE_vs_PDC_results.(xlsx|csv)` | `./figures/fig7/phylum_boxplots.(pdf|svg|tiff)` |

> Some figures in the manuscript (e.g., parts of Fig. 2b/5a/5b) were produced on **Origin** or **NovoMagic**; the underlying data files are deposited here, and the platform-specific figure generation is noted in file comments.

---

## Data scope, privacy, and rights

- All data here are **non-identifiable** and labeled with anonymized volunteer IDs (e.g., `V01`, `V02`).  
- Please **do not** attempt re-identification or linkage attacks.  
- You may reuse the data under **CC BY 4.0** (see `DATA_LICENSE.txt`) with proper attribution.  
- **Polite request**: besides citation, please leave a short note (issue or email) if you substantially reuse/extend the dataset—this helps us track impact and answer questions early.

---

## How to cite

- **Software (code)**: Cite this repository (see the “Cite this repository” button on GitHub, generated from `CITATION.cff`).  
- **Data**: Cite the repository and specific dataset/file where relevant.  
- **Paper**: Once accepted, we will add the **DOI and URL** here; please cite the paper for conceptual/methodological claims.

**Temporary citation (until the paper DOI is available):**

> Gao Han & Liu Quanjun (2025). *PDC Sampler Condensation Theory & Bioinformatics (ASE vs PDC)* (Version 1.0). GitHub repository. Available at: [URL to this repo]

---

## Contact

- First author: **Gao Han** — 220224443@seu.edu.cn  
- Corresponding author: **Liu Quanjun** — lqj@seu.edu.cn

We welcome questions, suggestions, and pull requests that improve clarity and reproducibility.



