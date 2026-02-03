# Metagenomic Profiling of Antimicrobial Resistance in Wastewater from Metropolitan Cities of India
---
## Contents
---
1. Introduction
2. Code
3. Data
---
## Introduction

The study analyses wastewater samples from March 2022 to March 2024 from 19 locations in four metropolitan cities of India. The samples were colleceted monthly once (447 samples) and profiled using shotgun metagenomics. 
The samples were sequenced on Illumina NovaSeq 6000 platform. Relative taxonomic and ARG abundances were computed directly from the reads. Subsequently, reads were assembled into contigs, which were then binned into metagenome assembled genomes (MAGs). These MAGs were further analyzed for taxonomic classification, ARG identification, profiling of MGEs and co-abundance network analysis.

## Code

This R script "main.R" can be used to generate images for the project "Metagenomic Profiling of Antimicrobial Resistance in Wastewater from Metropolitan Cities of India".
All the packages needed to run the code is listed at the top of the R script and should be installed before running the code.
The session information for the R enviroment used to generate the code is provided in "sessionInfo" file and should be used to get the package information.

## Data

All the data needed for the code is deposited at the [Zenodo](https://doi.org/10.5281/zenodo.17590646). 
The data should be downloaded in the same directory as the main.R file which should also be set as the working directory.
These data are in tabular format and contains processed data received from the specific tools.
Most of the tools give sample specific results, which were then read and combined in a tabular format. Any processing, including normalization is done within the code.

| file name | Description |
| ----------- | ----------- |
| k2_nt_20230502.inspect.species.krona.phylotree.txt | inspect file of k2_nt database downloaded on 2nd May 2023 |
| meta.tsv | Meta file for the data describing sampling city and collection month |
| bacterial_read_df.tsv | File with total read count (kraken) and domain level read count (bracken)|
| bracken_phylum.tsv | Bracken estimated read at the phyla level |
| bracken_genus.tsv | Bracken estimated read at the genra level |
| bracken_species.tsv | Bracken estimated read at the species level |
| rgi_bowtie.tsv | RGI output file |
| checkm2_quality_report.csv | CheckM2 quality report for all the MAGs |
| genomad_plasmid_HQ_MAG.tsv | geNomad output containing plasmid score for high quality MAGs |
| bat_HQ_MAG.tsv | BAT output for high quality MAGs classification |
| amr_HQ_MAG.tsv | BacAnt output for hiqh quality MAGs identifying AMR signatures |
| transposon_HQ_MAG.tsv | BacAnt output for hiqh quality MAGs identifying transposon signatures |
| sample_to_dereplicated_bin_coverage.tsv | coverage of HQ MAG for each samples |



