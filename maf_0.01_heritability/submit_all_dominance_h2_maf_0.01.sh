#!/bin/bash

# Ensure that the ldscore files and M_5_50, M have been copied up to the cluster.
# scp ~/Repositories/ldscgxe/1000G_EUR_Phase3_plink/ldscores_maf_0.01/* dpalmer@login01:/stanley/genetics_storage/analysis/ukbb_dominance/ldscores_ldscgxe/ldscores_maf_0.01_nov_2022/

# Run it
sumstats_folder=/psych/genetics_data/dpalmer/UKbb/ldsc-dominance-export
grep=$1

path_to_output="/psych/genetics_data/dpalmer/UKbb/h2_maf_0.01_results"
mkdir $path_to_output;

for sumstat_filename in $sumstats_folder/${grep}*both_sexes*; do
	phenotype=`echo $sumstat_filename | sed -e "s/^[\.\.//]*//" | sed -e "s/\..*//" | sed -e "s/^.*\///"`;
	echo "$phenotype submitted to the cluster for dominance heritability estimation";
	qsub -v sumstat_file="$sumstat_filename",phenotype="$phenotype",path_to_out="${path_to_output}" h2_dominance_template_maf_0.01.sh
done
