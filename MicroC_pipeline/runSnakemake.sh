#! /bin/bash -l
 
#SBATCH --partition=scu-cpu   # cluster-specific
#SBATCH --time=72:00:00   # HH/MM/SS
#SBATCH --mem=256Mb   # memory requested, units available: K,M,G,T
#SBATCH --mail-user=hpc_alert_ukl4001@outlook.com
#SBATCH --mail-type=ALL
#SBATCH -e snakeSlurm.err.%J
#SBATCH -o snakeSlurm.out.%J

source ~/.bashrc
cd $SLURM_SUBMIT_DIR
echo "SLURM_SUBMIT_DIR = $SLURM_SUBMIT_DIR"
echo "PWD at runtime    = $(pwd)"

###### Activating conda environment for snakemake
###### snakemake package should be installed in this conda enviromnment
conda activate snakemake

###### Receiving snakeFile from command line input
snakeFile=$1


###### Creating folders for output
if [[ ! -d ./log/ ]]; then
	mkdir ./log/
fi

if [[ ! -d ./log/Err_Out/ ]]; then
	mkdir ./log/Err_Out/
fi

##### Checking if snakefile is profived as input
if [[ ! -f ${snakeFile} ]];
then
	echo "snakefile does not exist!"
	exit 1
fi


##### Running snakefile 
echo "Using Snakefile: ${snakeFile}"
mkdir -p /athena/cayuga_0007/scratch/ukl4001/temp
export TMPDIR=/athena/cayuga_0007/scratch/ukl4001/temp

snakemake -s ${snakeFile} \
	--rerun-incomplete \
	--cluster-config config/slurmConfig.json \
	--use-conda \
	--conda-frontend conda \
	--latency-wait 200 \
	--printshellcmds \
	--cluster  \
	"sbatch --partition={cluster.partition} -J {rule} -o log/Err_Out/slurm-%j.out -e log/Err_Out/slurm-%j.err -N1 --cpus-per-task {cluster.threads} --time {cluster.time} --mem={cluster.mem}" \
	--jobs 500 \


echo "JOB FINISHED"
exit
