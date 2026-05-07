# WGBS Analysis - Homework 2

## Project Description
Whole Genome Bisulfite Sequencing (WGBS) analysis of esophageal squamous cell carcinoma (ESCC) tumour and matched normal tissue samples to identify differentially methylated cytosines (DMCs) and regions (DMRs).

## Samples
| Sample ID | Condition | GSM ID |
|-----------|-----------|--------|
| SRR11647648 | Normal | GSM4505856 |
| SRR11647649 | Normal | GSM4505857 |
| SRR11647658 | Tumor | GSM4505866 |
| SRR11647659 | Tumor | GSM4505867 |

## Data Description
- Paired-end WGBS data
- Sequenced on Illumina HiSeq2500
- Reference genome: GRCh37/hg19

## Software and Packages Used
- FastQC v0.12.1
- MultiQC v1.33
- Trim Galore v0.6.10
- BSMAP v2.89
- samtools
- MethylDackel
- R packages: ggplot2, corrplot, data.table, patchwork

## Analysis Steps
1. Raw data QC (FastQC, MultiQC)
2. Adapter trimming (Trim Galore)
3. Read mapping to GRCh37 (BSMAP)
4. Duplicate marking (samtools markdup)
5. Methylation bias assessment (MethylDackel mbias)
6. Methylation extraction on chromosome 20 (MethylDackel extract)
7. Coverage QC and statistical analysis (R)

## Reproducibility
All analysis was performed on a remote Linux server. Raw data downloaded from ENA using wget.
