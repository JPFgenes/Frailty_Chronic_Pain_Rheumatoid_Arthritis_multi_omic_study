# Frailty, chronic pain and rheumatoid arthritis: analysis methods

This repository provides the following scripts: bidirectional two-sample Mendelian randomisation (MR), pairwise GWAS (GWAS-PW), shared-variant credible-set summaries, colocalisation, cross-sectional polygenic-score (PRS)/protein-EpiScore associations and networks, and longitudinal protein–frailty models.

The scripts are intended to **illustrate and allow adaptation of the statistical workflow using appropriately prepared data**. They are not an end-to-end reproduction of the published numerical estimates. No participant-level data are included.

## Contents

| File | Purpose |
| --- | --- |
| `01_pw_gwas_pipeline.sh` | Run GWAS-PW, extract variant posteriors, call postprocessing and optionally annotate nearest genes using BEDtools. |
| `02_pw_gwas_postprocess.R` | Summarise regional Model 0–4 probabilities, select regions by PPA thresholds and construct nominal 90% shared-variant credible sets. |
| `03_colocalisation.R` | Perform `coloc.abf` in selected PW-GWAS regions, with allele alignment, variant QC and quantitative/case-control settings. |
| `04_EpiScore_PRS_Network.R` | Standardise PRSs; run cross-sectional PRS–biomarker and EpiScore–FI regressions; apply FDR; plot PRS–EpiScore–FI and protein correlation networks. |
| `05_EpiScore_Longitudinal.R` | Fit FI and EpiScore/protein mixed-effects trajectories; extract participant slopes; test protein slope associations with frailty progression. |
| `06_Bidirectional_MR.R` | Select/clump instruments; harmonise and Steiger-filter; estimate bidirectional MR effects and sensitivity results; apply IVW FDR. |

## Software

- R (manuscript: 4.3.1), with `data.table`, `coloc` (manuscript: 5.2.3), `TwoSampleMR`, `igraph` and `lme4`; `MRPRESSO` is needed for MR-PRESSO sensitivity analysis.
- GWAS-PW (the pipeline expects an installation supporting `--extract-snp-posteriors`), PLINK 1.9 and a matching European-ancestry LD reference panel for MR; BEDtools is optional for gene annotation.
- GRCh37/hg19-aligned inputs for the GWAS-based analyses; appropriate LD blocks for GWAS-PW.

Install tools and dependencies separately. Paths in scripts are placeholders; R packages, binary versions and inputs must be chosen for the research dataset.

## Required upstream data preparation

The scripts assume that the following inputs have already been generated and quality-controlled, using the manuscript and cited sources for details:

- **Frailty Index (FI):** cohort-specific health deficits have been coded and aggregated as an FI.
- **Protein EpiScores:** prediction models previously trained on DNA methylation have been projected to participants; either supply technically adjusted EpiScores, or explicitly enable and review the example residualisation in script `04`.
- **PRSs:** scores have been calculated upstream (the study used SBayesR); script `04` standardises supplied PRSs within each cohort/wave.
- **GWAS data:** summary statistics have been harmonised to GRCh37, appropriately allele oriented, and reshaped into each script's required columns. GWAS-PW input preparation, including the necessary Z-statistic/sample-overlap setup, is upstream of script `01`.

## 1. GWAS-PW, credible sets and annotation

Run from the repository root after installing `gwas-pw` and supplying an appropriate two-trait GWAS-PW input and an LD-block BED:

```bash
bash 01_pw_gwas_pipeline.sh \
  data/gwaspw_pair.tsv reference/ld_blocks.bed TraitA TraitB \
  results/traitA_traitB NONE FALSE
```

The sixth argument may be a four-column GRCh37 gene BED (`chr`, zero-based `start`, half-open `end`, `gene_name`) instead of `NONE`; it enables nearest-gene assignment via BEDtools. The last argument toggles an **additional** APOE exclusion (`TRUE` or `FALSE`). The pipeline excludes extended MHC blocks (chr6:25–35 Mb) during postprocessing.

The postprocessor uses PPA3 >= 0.50 for suggestive shared regions, PPA3 >= 0.80 for nominal credible sets, and produces PPA3 >= 0.90 and PPA4 >= 0.90 block lists for downstream coloc. Its CS90 algorithm normalises extracted variant-level `PPA_3` weights *within each block* and includes the SNP reaching 90% cumulative weight. The interpretation of these extracted values as calibrated posterior weights must be checked against the installed GWAS-PW implementation. Some original analyses used a different crossing-SNP rule, and this change can alter set size. BEDtools is used for **annotation**, not for the Bayesian inference itself.

## 2. Colocalisation

**Unlike the other R scripts, `03_colocalisation.R` requires a one-row, tab-separated configuration file.** This file is not included in the six-script package. To run coloc, create `config/coloc_config.tsv` with these columns:

```text
trait1_file trait2_file blocks_file trait1_type trait2_type trait1_sdY trait2_sdY trait1_N trait2_N trait1_case_fraction trait2_case_fraction output_prefix p12
```

Use **tabs**, not spaces, between the fields and add one data row containing your own paths and trait metadata. Set `trait*_type` to `quant` or `cc`. For quantitative traits, provide a defensible `sdY`, or provide `N` and observed MAF for the script to estimate it; for case-control traits, provide total sample size (`N`) and **sample case fraction** (the proportion of cases in the analysed GWAS). Specify `p12` explicitly and ensure the choice is consistent with the results being compared. The file `blocks_file` must point to the `_blocks_filtered.tsv` generated in step 1 (columns `block,chr,st,sp,P3,P4`). GWAS inputs require `chr,pos,snp,EA,OA,beta,se,p,MAF`.

```bash
Rscript 03_colocalisation.R config/coloc_config.tsv
```

The script selects blocks with PPA3 >= 0.90 or PPA4 >= 0.90, expands boundaries by 250 kb, aligns directly matching/swapped alleles, removes ambiguous A/T and C/G SNPs, and requires observed MAF >= 0.01 in both traits, p < 0.05 in trait 1 and at least 100 QC SNPs. It does not handle every strand-complement case; inputs should be pre-harmonised. The original analyses included different `p12` settings in exploratory runs, so check the prior against final manuscript tables.

## 3. Cross-sectional EpiScore, PRS and network analyses

Edit the **USER SETTINGS** block near the top of `04_EpiScore_PRS_Network.R` (input path; participant/cohort/wave fields; FI, PRS, EpiScore and optional measured-protein columns). Its input table needs one row per participant/wave, age, sex and four genetic PCs. It fits separate PRS–biomarker linear models adjusting for age, sex and PCs and EpiScore–FI models adjusting for age and sex (minimum 50 complete observations). It applies BH FDR and creates directed PRS–EpiScore–FI and optional protein–protein Pearson correlation networks (|r| > 0.30).

```bash
Rscript 04_EpiScore_PRS_Network.R
```

## 4. Longitudinal EpiScore/protein–frailty analyses

Edit the **USER SETTINGS** block of `05_EpiScore_Longitudinal.R` and supply one row per participant and wave. The script models protein values (up to wave 4) and FI (up to wave 5) with random intercepts and age slopes:

```r
protein ~ age + sex + (1 + age | ID)
FI      ~ age + sex + (1 + age | ID)
```

It extracts participant-level slopes, excludes biomarker slopes with between-person SD < 0.01, standardises protein slopes, then fits:

```r
FI_slope ~ protein_slope_z + protein_intercept + baseline_age + sex
```

and BH-adjusts the protein association p-values.

```bash
Rscript 05_EpiScore_Longitudinal.R
```

The longitudinal input and age units must reflect the intended study design. In particular, verify age centring, empirical Bayes intercept interpretation and slope scaling against the original analysis before comparing exact effect sizes.

## 5. Bidirectional Mendelian randomisation

Edit the **USER SETTINGS** block at the top of `06_Bidirectional_MR.R`. Provide one GWAS table per trait with `SNP,chr,pos,ea,oa,beta,se,p,eaf,n`, where `eaf` is **effect-allele frequency**, not minor-allele frequency. Binary traits additionally require `ncase,ncontrol`, and a justified **population prevalence** for Steiger directionality (distinct from the GWAS sample case fraction used in coloc). Specify forward/reverse trait pairs, PLINK 1.9 binary and a compatible 1000 Genomes EUR PLINK reference.

```bash
Rscript 06_Bidirectional_MR.R
```

The template uses p < 5e-8, local LD clumping (r² < 0.001, 10 Mb), harmonisation (`action=2`), F-statistic > 10, Steiger filtering, IVW and MR-Egger/median/mode sensitivity methods, heterogeneity, Egger intercept, leave-one-out, optional MR-PRESSO, and BH FDR for IVW estimates. **The current script defaults to `exclude_mhc <- TRUE`**. All trait metadata must be valid for Steiger filtering.
