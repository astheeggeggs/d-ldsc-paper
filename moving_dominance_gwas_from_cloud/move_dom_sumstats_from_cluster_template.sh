#\!/bin/bash
#$ -N get_ldsc_file
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=4g

source /broad/software/scripts/useuse
use R-3.5
export PATH="/home/unix/dpalmer/google-cloud-sdk/bin/:$PATH"

echo "gs_file = ${gs_file}"
gsutil -m cp ${gs_file} /stanley/genetics_storage/analysis/ukbb_dominance/sumstats-dominance-export/both_sexes/
