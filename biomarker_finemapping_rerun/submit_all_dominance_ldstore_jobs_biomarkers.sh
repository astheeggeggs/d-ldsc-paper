#!/bin/bash 
# Run LD-store jobs (using dsub) from the cluster

maf_locus=0.05
maf_finemap=0
p_hwe_finemap=1e-10

# biomarkers
dsub_submission_files=/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/biomarkers/inputs/z_and_task_files_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}/dominance

cd ${dsub_submission_files}

for job in ${dsub_submission_files}/*ldstore.dsub*; do
	echo "${job}";
	bash "${job}";
done
