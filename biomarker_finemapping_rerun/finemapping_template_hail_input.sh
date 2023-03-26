#\!/bin/bash
#$ -N get_y
#$ -cwd
#$ -o /psych/genetics_data/dpalmer/logs/
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=32g
#$ -l h_rt=1:00:00
#$ -l os=RedHat7

source /broad/software/scripts/useuse

export PATH="/psych/genetics_data/dpalmer/anaconda2/bin:$PATH"
export PATH="/home/unix/dpalmer/google-cloud-sdk/bin/:$PATH"

# conda env create --file environment.yml
# conda env create --file finemapping_input_environment.yml
source activate finemap_inputs
use R-3.5

# Removed sex as an option as we're just initially focusing on the both_sexes GWASes.
bash /home/unix/dpalmer/Repositories/ldscgxe_dominance_paper/biomarker_finemapping_rerun/finemapping_hail_input.sh ${pheno} ${full_phenotype_file} ${incl_file} ${dom_file} ${add_file} ${out} ${maf_locus} ${maf_finemap} ${p_hwe_finemap}
