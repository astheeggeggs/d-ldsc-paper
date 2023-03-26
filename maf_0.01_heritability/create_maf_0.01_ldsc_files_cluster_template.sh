#\!/bin/bash
#$ -N get_ldsc_file
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=8g

source /broad/software/scripts/useuse
use R-3.5
export PATH="/home/unix/dpalmer/google-cloud-sdk/bin/:$PATH"

echo "batch = ${batch}"
Rscript --vanilla  create_maf_0.01_ldsc_files_cluster.r ${batch}
