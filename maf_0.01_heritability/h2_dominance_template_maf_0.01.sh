#\!/bin/bash
#$ -N run_h2_pheno
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=8g

source /broad/software/scripts/useuse
use GCC-5.2
use GSL
use UGER
use Anaconda
reuse -q PLINK
reuse -q .zlib-1.2.8
reuse -q PSEQ

source activate d-ldsc

ldscore_files="/stanley/genetics_storage/analysis/ukbb_dominance/ldscores_ldscgxe/ldscores_maf_0.01_nov_2022/1000G.EUR.QC.chr@";
ldscore_files_weights="/stanley/genetics_storage/analysis/ukbb_dominance/ldscores_ldscgxe/ldscores_maf_0.01_nov_2022/1000G.EUR.QC.chr@";
out="${path_to_out}/add_dom_maf_0.01.${phenotype}" 
echo "./get_h2.py --dominance --ref-ld-chr $ldscore_files --w-ld-chr $ldscore_files_weights --n-blocks 200 --write-h2 --h2 $sumstat_file --chisq-max 10000000000 --out $out --pheno-name $phenotype"
python /home/unix/dpalmer/Repositories/d-ldsc/get_h2.py --additive --dominance --ref-ld-chr $ldscore_files --w-ld-chr $ldscore_files_weights --n-blocks 200 --write-h2 --h2 $sumstat_file --chisq-max 10000000000 --out $out --pheno-name $phenotype
