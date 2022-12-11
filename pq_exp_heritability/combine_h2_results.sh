#!/bin/bash

path_to_h2_folder="/psych/genetics_data/dpalmer/h2-pq-exp"
outdir=${path_to_h2_folder}/combined
mkdir $outdir

for type in {"add","dom"}; do
	for alpha in {"-1.0","-0.5","-0.25","0.0","0.25","0.5","1.0"}; do
		file_to_write="${outdir}/${type}.alpha.${alpha}.h2";
		echo $file_to_write;
		echo -n > $file_to_write;
		for h2_filename in "${path_to_h2_folder}/${type}.alpha.${alpha}.*h2"; do
			# phenotype=`echo $h2_filename | sed -e "s/^[\.\.//]*//" | sed -e "s/\..*//" | sed -e "s/^.*\///"`;
			echo $h2_filename
			cat $h2_filename >> $file_to_write;
		done
		sed -i '1iphenotype\tmean_chi2\tlambdaGC\tintercept\tintercept_se\tintercept_z\tintercept_p\tratio\tratio_se\th2_observed\th2_observed_se\th2_liability\th2_liability_se\th2_z\th2_p' $file_to_write
	done
done
