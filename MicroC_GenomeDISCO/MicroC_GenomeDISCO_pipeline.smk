###############################################################################
### Running GenomeDISCO                                                     ###
###############################################################################

### Importing python packages
import os
import pandas as pd
import re

### Configuration
configfile: "config/config.json"

###############################################################################
# RULES
rule all:
	input:
		expand("/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/log/genomedisco_{num}.ok", num = list(range(58, 63)))

###############################################################################

rule genomedisco:
	conda: "genomedisco_python2"
	threads: 4
	input:
		samples = "/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/metadata/metadata_{num}.samples",
		pairs = "/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/metadata/metadata_{num}.pairs",
	output:
		success = "/athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/log/genomedisco_{num}.ok"
	params:
		num = "{num}",
		outdir = "/athena/apostoloulab/scratch/ukl4001/data/GenomeDSICO/comparison_{num}"
	shell:
		"""
		mkdir -p {params.outdir}

		genomedisco run_all --metadata_samples {input.samples} \
					 --metadata_pairs {input.pairs} \
					 --bins /athena/apostoloulab/scratch/ukl4001/reference/mm10/mm10.bin.50kb.bed.gz \
					 --parameters_file /athena/apostoloulab/scratch/ukl4001/MtoG1_analysis_code/MicroC_GenomeDISCO/parameters.txt \
					 --outdir {params.outdir}

		# Made it?
		touch {output.success}
		"""
