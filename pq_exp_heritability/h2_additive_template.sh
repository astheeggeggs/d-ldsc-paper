#\!/bin/bash
#$ -N run_h2_pheno
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=4g

source /broad/software/scripts/useuse
use GCC-5.2
use GSL
use UGER
use Anaconda
reuse -q PLINK
reuse -q .zlib-1.2.8
reuse -q PSEQ

source activate d-ldsc

for alpha in {"-1.0","-0.5","-0.25","0.0","0.25","0.5","1.0"}; do
	ldscore_files="/psych/genetics_data/dpalmer/ldscores-pq-exp/ldscores/1000G.EUR.QC.alpha.${alpha}.chr@";
	ldscore_files_weights="/psych/genetics_data/dpalmer/ldscores-pq-exp/ldscores/1000G.EUR.QC.alpha.${alpha}.chr@";
	out="${path_to_out}/add.alpha.${alpha}.${phenotype}" 
	echo "./get_h2.py --additive-orig --ref-ld-chr $ldscore_files --w-ld-chr $ldscore_files_weights --n-blocks 200 --write-h2 --h2 $sumstat_file --chisq-max 10000000000 --out $out --pheno-name $phenotype"
	python /home/unix/dpalmer/Repositories/d-ldsc/get_h2.py --additive-orig --ref-ld-chr $ldscore_files --w-ld-chr $ldscore_files_weights --n-blocks 200 --write-h2 --h2 $sumstat_file --chisq-max 10000000000 --out $out --pheno-name $phenotype
done
