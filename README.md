# Post-mitotic transcriptional activation and 3D regulatory interactions show locus- and differentiation-specific sensitivity to cohesin depletion

This repository contains the analysis code for the following publication:

> **Post-mitotic transcriptional activation and 3D regulatory interactions show locus- and differentiation-specific sensitivity to cohesin depletion**
>
> UkJin Lee, Alejandra Laguillo-Diego, Daniela Magliulo, Wilfred Wong, Kritika Kasliwal, Zhangli Ni, Lingling Cheng, Jieru Li, Bobbie Pelham-Webb, Alexandros Pertsinidis, Christina Leslie, Effie Apostolou
>
> Preprint: [bioRxiv 2025.02.13.638153](https://www.biorxiv.org/content/10.1101/2025.02.13.638153v1)

## Overview

This study investigates how cohesin depletion affects transcriptional reactivation and 3D genome organization during the mitosis-to-G1 transition in mouse embryonic stem cells. Analyses span RNA-seq, ChIP-seq, Micro-C, 4C-seq, and integration with published datasets.

## Repository structure

### Micro-C processing pipeline

| Directory | Description |
|-----------|-------------|
| `MicroC_pipeline/` | 5-step Snakemake pipeline: FASTQ QC and BWA-MEM alignment, pairtools dedup, .hic/.mcool generation, KR normalization, and expected contact calculation |
| `MicroC_downsampling/` | Downsample all samples to equal read depth using cooltools |
| `MicroC_GenomeDISCO/` | Reproducibility scoring between replicates and conditions using GenomeDISCO |

### Micro-C analysis

| Directory | Description |
|-----------|-------------|
| `MicroC_compartment/` | A/B compartment calling using CscoreTool |
| `MicroC_TAD/` | TAD boundary calling using cooltools insulation score |
| `MicroC_CvD/` | Contact-vs-distance P(s) decay curves |
| `MicroC_microcompartment/` | Microcompartment detection using Chromosight |

### Loop analysis

| Directory | Description |
|-----------|-------------|
| `MicroC_loop_chromosight/` | Loop detection using Chromosight |
| `MicroC_loop_fithic/` | Loop detection using FitHiC |
| `MicroC_loop_hicdcp/` | Loop detection using HiCDC+ |
| `MicroC_loop_analysis/` | Loop comparison across callers, APA pileups (coolpuppy), GO/LOLA enrichment, and other miscellaneous loop analysis|
| `MicroC_loop_Hansen/` | Re-scoring loops from [Jusuf et al.](https://doi.org/10.1101/2025.01.13.632736) on internal Micro-C data |
| `MicroC_loop_Tjian/` | Re-scoring loops from [Hsieh et al.](https://doi.org/10.1038/s41588-022-01223-8) on internal Micro-C data |

### Virtual 4C and 4C-seq

| Directory | Description |
|-----------|-------------|
| `MicroC_v4C/` | Virtual 4C extraction from Micro-C .pairs files with sliding-window coverage |
| `4C_postprocessing/` | Post-processing of conventional 4C-seq data |
| `4C_v4C_visualization/` | Combined visualization of 4C and virtual-4C tracks |

### RNA-seq analysis

| Directory | Description |
|-----------|-------------|
| `RNA_RAD21/` | Integration of internal G1-phase RAD21 depletion data with published async |
| `RNA_NIPBL/` | Analysis of published NIPBL depletion time-course RNA-seq |

### Figures and utilities

| Directory | Description |
|-----------|-------------|
| `Miscellaneous_figure/` | Additional figure panels and ChIP-seq track visualization |
| `script/` | Shared R functions (RNA-seq analysis, figure rendering, parameters), shell helpers, and Python utilities |
| `reference_genome/` | mm10 gene annotation generation from Gencode vM25 via biomaRt |
| `figure/` | Output directory for generated figures |
| `data/` | Placeholder for input data files |
| `r_data/` | Placeholder for intermediate R data objects |

## Reference genome

All analyses use the **mm10 (GRCm38)** mouse genome assembly with **Gencode vM25** gene annotations.

## Data availability

Raw and processed data are available on the GEO database under accession number **[GSE288339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE288339)**.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
