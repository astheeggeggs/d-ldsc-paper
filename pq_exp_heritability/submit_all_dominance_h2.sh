#!/bin/bash

sumstats_folder=$1; 		# /stanley/genetics_storage/analysis/ukbb_dominance/ldsc-dominance-export/sumstats-files-1, 2, etc;
path_to_output=$2; 			# /psych/genetics_data/dpalmer/h2-pq-exp/;

mkdir $path_to_output;

for sumstat_filename in $sumstats_folder/*both_sexes*; do
	phenotype=`echo $sumstat_filename | sed -e "s/^[\.\.//]*//" | sed -e "s/\..*//" | sed -e "s/^.*\///"`;
	echo "$phenotype submitted to the cluster for dominance heritability estimation";
	qsub -v sumstat_file="$sumstat_filename",phenotype="$phenotype",path_to_out="${path_to_output}" h2_dominance_template.sh
done
