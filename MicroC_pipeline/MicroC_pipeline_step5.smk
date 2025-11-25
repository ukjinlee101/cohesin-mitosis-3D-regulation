###############################################################################
### Micro-C pipeline for PE dataset                                         ###
###############################################################################

###############################################################################

# REFERENCE
# https://micro-c.readthedocs.io/en/latest/before_you_begin.html
# .hic and .mcool file have KR, SCALE, VC, VC_SQRT normalization column. However, .mcool will not have 'sum' metadata calculated.
# hicConvertFormat extracts cool file of specific resolution from mcool file and calculate weight column based on the KR column. Weight = 1/KR.
# However, for consistency with other dataset, KR column will be calculated again using cooler balance. KR balancing for Jucier and Cooler are therotically same, but numerically different.
###############################################################################


### Importing python packages
import os
import pandas as pd
import re

### Configuration
configfile: "config/config.json"
sampleBatchList = ["G1DMSO_B1", "G1DMSO_B2",
				   "G1dTAG_B1", "G1dTAG_B3",
				   "G1A485_B2", "G1A485_B3"]
resList = [25000, 10000, 5000]
###############################################################################
# RULES
rule all:
	input:
		expand("result/log/success/normCoolBatch_{sampleBatch}_{res}.ok", sampleBatch = sampleBatchList, res = resList),


###############################################################################
rule normCoolBatch:
	conda: "hicexplorer"
	threads: 16
	input:
		mcool = "/athena/apostoloulab/scratch/ukl4001/data/mcool_batch/{sampleBatch}_allRes.mcool"
	output:
		success = "result/log/success/normCoolBatch_{sampleBatch}_{res}.ok",
		cool = "/athena/apostoloulab/scratch/ukl4001/data/cool_norm_batch/{sampleBatch}_{res}bp_KR.cool",
		exp = "/athena/apostoloulab/scratch/ukl4001/data/cool_norm_batch/{sampleBatch}_{res}bp_KR_exp.tsv"
	params:
		res = "{res}"
	shell:
		"""

		hicConvertFormat -m {input.mcool}:://resolutions/{params.res} --inputFormat cool \
			--outputFormat cool -o {output.cool} --correction_name KR

		cooler balance --force -p 16 {output.cool}

		cooltools expected-cis -p 16 -o {output.exp} {output.cool}

		# Made it?
		touch {output.success}
		"""

