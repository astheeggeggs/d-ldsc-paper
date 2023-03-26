#\!/bin/bash
#$ -N run_h2_sim
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=4g

source /broad/software/scripts/useuse
# use GCC-5.2
# use GSL
# use UGER
# use Anaconda
# reuse -q PLINK
# reuse -q .zlib-1.2.8
# reuse -q PSEQ

source activate d-ldsc

ldscore_files="/stanley/genetics_storage/analysis/ukbb_dominance/ldscores_ldscgxe/1000G_EUR_Phase3_hm3_chr@"
ldscore_files_weights="/stanley/genetics_storage/analysis/ukbb_dominance/ldscores_ldscgxe/1000G_EUR_Phase3_hm3_chr@"
sumstats_folder="/psych/genetics_data/dpalmer/UKbb/ldsc-dominance-biomarkers-export-rerun/"
out="/psych/genetics_data/dpalmer/UKbb/ldsc-dominance-biomarkers-export-rerun/h2-results"

mkdir $out;

for sumstat_filename in $sumstats_folder/*; do
	phenotype=`echo $sumstat_filename | sed -e "s/^[\.\.//]*//" | sed -e "s/\..*//" | sed -e "s/^.*\///"`;
	echo "$phenotype submitted to the cluster for heritability estimation";
	python /home/unix/dpalmer/Repositories/d-ldsc/get_h2.py --additive --dominance --ref-ld-chr $ldscore_files --w-ld-chr $ldscore_files_weights --n-blocks 200 --write-h2 --h2 $sumstat_filename --chisq-max 10000000000 --out $out/$phenotype --pheno-name $phenotype
done
