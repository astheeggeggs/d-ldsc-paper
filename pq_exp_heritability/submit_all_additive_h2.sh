#!/bin/bash

sumstats_folder=$1; 		# /stanley/genetics/analysis/ukbb_sumstats/ldsc_sumstats_additive/round2
path_to_output=$2; 			# /psych/genetics_data/dpalmer/h2-pq-exp/;
grep=$3;					# What to grep for in the folder, because there's a massive number of files and don't want to submit too much.

mkdir $path_to_output;

for sumstat_filename in $sumstats_folder/${grep}*both_sexes*; do
	phenotype=`echo $sumstat_filename | sed -e "s/^[\.\.//]*//" | sed -e "s/\..*//" | sed -e "s/^.*\///"`;
	echo "$phenotype submitted to the cluster for additive heritability estimation";
	qsub -v sumstat_file="$sumstat_filename",phenotype="$phenotype",path_to_out="$path_to_output" h2_additive_template.sh
done
