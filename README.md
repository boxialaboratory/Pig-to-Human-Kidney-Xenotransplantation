# Pig-to-Human Kidney Xenotransplantation: scRNA-seq Data Analysis

This repository contains the R scripts used for analyzing scRNA-seq, longitudinal bulk RNA-seq, and longitudinal host peripheral blood mononuclear cells sequencing data derived from pig-to-human kidney xenotransplantation conducted in 2021. Each script is focused on a different aspect of the data analysis, providing a comprehensive look into the xenotransplantation process and its outcomes.

## Contents
1. [Overview](#overview)
2. [Associated Publication](#associated-publication)
3. [Scripts](#scripts)
4. [How to Run the Analysis](#how-to-run-the-analysis)

## Overview
Xenotransplantation is a promising solution to address the challenge of organ donor shortage. Two porcine-to-human kidney xenotransplantations were performed in 2021, yet the cellular physiology between xenograft and the recipient remains largely unknown. We conducted single-cell RNA sequencing on xenografts and their contralateral control kidneys, longitudinal bulk RNA-seq on xenograft core biopsies, and longitudinal single-cell RNA sequencing on peripheral blood mononuclear cells to dissect the time-resolved xenograft-recipient physiological interactions across time. Longitudinal and single-cell RNA-seq analyses of porcine kidneys and recipient’s PBMCs revealed time-resolved cellular dynamics of xenograft-recipient interactions. All data and scripts used in the analysis can be found in this respository.

## Associated Publication
The analyses contained in this repository are integral to our associated publication titled Cellular dynamics of pig-to-human kidney xenotransplantation, published in Med, 2024. In the paper, we detail the comprehensive single-cell RNA sequencing (scRNA-seq) and computational analyses conducted on the xenografts and their contralateral untransplanted kidneys, as well as on the peripheral blood mononuclear cells (PBMCs) of the recipients. These analyses allowed us to better understand the physiological impact of xenotransplantation, providing insights into cellular dynamics, immune responses and potential early signs of rejection.

## Scripts
The repository contains the following R scripts:

- `01_data_preprocessing.R`: Preprocessing and exploratory analysis of the scRNA-seq data .
- `02_immune_cell_analysis.R`: Analysis of the immune cells within the single cell RNA sequencing population and validation of gene expression trajectory with longitudinal bulk-RNA-seq data.
- `03_rejection_type_analysis.R`: Analysis xenografts rejections and validation of gene expression trajectory with longitudinal bulk-RNA-seq data.
- `04_proliferating_population_analysis.R`: Analysis of the proliferating proximal tubule cells in the xenograft.
- `05_xenograft_damage_analysis.R`: Differentially Expressed Genes analysis revealing global damage and cell type-specific repair of the xenograft.
- `06_pbmc_pre_processing.R`: Handles preprocessing of PBMC scRNA-seq data.
- `07_annotations_EDA.R`: Exploratory data analysis of of immune cell population in PBMC.
- `08_gene_clustering_analysis.R`: Performs time-course coexpression analysis in each immune cell population.
- `09_gene_set_analysis.R`: Analysis of co-expressed gene sets in each immune cell population to identify immune responses.

Each script is a standalone R file that performs a specific analysis and can be executed independently. Comments within each file provide necessary context and understanding of the process.

## How to Run the Analysis
Each script can be run independently in an R environment. Ensure that you have all the necessary R packages installed before running the scripts. Follow the numbered order of the scripts to replicate our analysis workflow.

This project is open for community access, contributions and collaborations are welcome. Please maintain the naming convention for consistency and easy understanding.
