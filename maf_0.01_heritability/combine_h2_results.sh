#!/bin/bash

path_to_h2_folder="/psych/genetics_data/dpalmer/UKbb/h2_maf_0.01_results"
outdir=${path_to_h2_folder}/combined
mkdir $outdir

file_to_write="${outdir}/add_dom_maf_0.01.h2";
echo $file_to_write;
echo -n > $file_to_write;
for h2_filename in "${path_to_h2_folder}/*h2"; do
	echo ${h2_filename};
	cat $h2_filename >> $file_to_write;
done
sed -i '1iphenotype\tmean_chi2_A\tlambdaGC_A\tintercept_A\tintercept_se_A\tintercept_z_A\tintercept_p_A\tratio_A\tratio_se_A\th2_observed_A\th2_observed_se_A\th2_liability_A\th2_liability_se_A\th2_z_A\th2_p_A\tmean_chi2_D\tlambdaGC_D\tintercept_D\tintercept_se_D\tintercept_z_D\tintercept_p_D\tratio_D\tratio_se_D\th2_observed_D\th2_observed_se_D\th2_liability_D\th2_liability_se_D\th2_z_D\th2_p_D' $file_to_write

# Moving from the cluster to local dir
scp dpalmer@login01:/psych/genetics_data/dpalmer/UKbb/h2_maf_0.01_results/combined/* h2_estimates/
