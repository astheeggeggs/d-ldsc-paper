#\!/bin/bash
#$ -N get_ldsc_file
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=8g

source /broad/software/scripts/useuse
# use GCC-5.2
# use GSL
# use UGER
# use Anaconda
# reuse -q PLINK
# reuse -q .zlib-1.2.8
# reuse -q PSEQ
use R-3.5

echo "batch = ${batch}"
Rscript --vanilla  create_maf_0.01_ldsc_files_cluster.r ${batch}
