# Supplementary code accompanying manuscript: A longitudinal single-cell atlas to predict outcome and toxicity after BCMA-directed CAR T cell therapy in multiple myeloma

This repository contains the code used to produce the results of Michael Rade, David Fandrei, Markus Kreuz et al. A longitudinal single-cell atlas to predict outcome and toxicity after BCMA-directed CAR T cell therapy in multiple myeloma, DOI XXXXXXXXXXXX

## Singularity

Most scripts for the single-cell analysis were developed in a Singularity image with Rstudio server. [See README](singularity/) in `./singularity/`. 

## Zenodo repository
Restricted access due to GDPR regulations and data can be accessed by placing a request via Zenodo (https://zenodo.org/records/14732727).

Includes:
- Singularity images for R/RStudio server and Souporcell
- All R packages used for publication
- Seurat object used for publication
- Cell ranger output

## Reproduction

``` sh
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Pre-Processing
# FASTQ -> Cell Ranger -> Souporcell -> Seurat
# The script must be executed in the base path (./) of this repo
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

### Sequencing run 1:
# see https://github.com/fraunhofer-izi/Rade_Grieb_et_al_2023
Rscript code/01_cohorts/ukl_batch_1/01_cellranger_to_seurat.R

### Sequencing run 2:
Rscript code/01_cohorts/ukl_batch_2/01_cellranger.R
Rscript code/01_cohorts/ukl_batch_2/02_souporcell.R
Rscript code/01_cohorts/ukl_batch_2/03_cellranger_to_seurat.R

### Sequencing run 3:
# The FASTQ files do not correspond to the naming convention of IlIumina. So we 
# have to rename them first.
Rscript code/01_cohorts/ukl_batch_3/01_fastq_symlinks.R

# In some cases, Ide and Cilta samples have been multiplexed. This means that
# you have to run cellranger twice, once with a custom Cilta reference and once
# with a custom Ide reference.
Rscript code/01_cohorts/ukl_batch_3/02_cellranger_ide.R

Rscript code/01_cohorts/ukl_batch_3/02_cellranger_cilta.R
Rscript code/01_cohorts/ukl_batch_3/03_souporcell.R
Rscript code/01_cohorts/ukl_batch_3/04_cellranger_to_seurat_ide.R
Rscript code/01_cohorts/ukl_batch_3/04_cellranger_to_seurat_cila.R

# Output is a Seurat object -> cells are now associated with donors
Rscript code/01_cohorts/ukl_batch_3/05_souporcell_demux.R

### Sequencing run 4:
Rscript code/01_cohorts/ukl_batch_4/03_cellranger_to_seurat.R
Rscript code/01_cohorts/ukl_batch_4/04_souporcell_demux.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Pre-Processing
# merge Objects -> QC -> Cell annotation -> Integration
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript code/01_cellranger.R
Rscript code/02_qc_merge_seurat.R
Rscript code/03_anno_1.R
Rscript code/04_anno_2.R
Rscript code/05_integration.R

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Main Figures
# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
Rscript publication/figure_scripts/fig_01.R
Rscript publication/figure_scripts/fig_02.R
Rscript publication/figure_scripts/fig_03_01_slingshot.R
Rscript publication/figure_scripts/fig_03_02.R
Rscript publication/figure_scripts/fig_04_02.R
Rscript publication/figure_scripts/fig_05.R
Rscript publication/figure_scripts/fig_06.R


```
