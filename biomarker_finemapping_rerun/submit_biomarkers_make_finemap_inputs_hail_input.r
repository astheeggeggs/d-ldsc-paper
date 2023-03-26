library(data.table)
library(dplyr)

maf_locus <- 0.05
maf_finemap <- 0
p_hwe_finemap <- 1e-10

# Go through each of the the different classes of phenotype, loop over them including the required phenotype and summary stat information

# First, check to ensure that all of the phenotype files on the cluster are the same as those on the cloud.

# Full phenotypes files
biomarkers_cluster <- "/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/ukb30163.restricted_and_phesant_formatted.biomarkers.tsv.gz"

incl_folders <- "/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/inputs/incl_files/biomarkers"
out_folders <- "/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/biomarkers"

awk_script <- "~/Repositories/ldscgxe/dominance_fine_mapping/incl_file_preparation.awk"

# phenotype_dt <- fread("gsutil cat gs://ukbb-dominance/phenotypes_with_significant_dominance_hit_no_cts_raw.tsv.bgz | zcat")
phenotype_dt <- fread(cmd = "gsutil cat gs://ukbb-dominance/phenotypes_with_significant_dominance_hit_no_cts_raw_0.05_pHWE_1e-6.tsv.bgz | zcat")
phenotypes <- phenotype_dt$phenotype
# Filter this list down to the 1,060 phenotypes in the curated list
# This file was created using generate_raymonds_primary_both_sex_phenotypes.r
# gsutil cp ~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv gs://ukbb-dominance/
phenotype_filter <- fread(cmd = "gsutil cat gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
phenotypes <- intersect(
	unlist(phenotype_filter %>% filter(dominance_qc_no_ordinal) %>% select(phenotype)),
	phenotypes
	)

for(i in 1:length(phenotypes)) {
	# Check to see if this phenotype is in this full phenotype file.
	file_phenotypes <- strsplit(system(paste("zcat", biomarkers_cluster, "| head -1"), intern=TRUE)[[1]], split='\t')[[1]]
	
	# Different behavior for the biomarkers (stored in a subfolder).
	file_phenotypes <- paste0(file_phenotypes, "_irnt")
	dom_file <- paste0("/psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/biomarkers_rerun/both_sexes/", phenotypes[i], ".gwas.imputed_v3.both_sexes.rerun.tsv.bgz")
	add_file <- paste0("/stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/biomarkers/", phenotypes[i], ".gwas.imputed_v3.both_sexes.tsv.bgz")

	if(phenotypes[i] %in% file_phenotypes) {
		cat("Submitting phenotype:", phenotypes[i], "\n")
		system(paste0('qsub -v ',
			'pheno=', phenotypes[i], ',',
			'full_phenotype_file=', biomarkers_cluster, ',',
			'incl_file=', incl_folders, '/', phenotypes[i], '.incl,',
			'dom_file=', dom_file, ',',
			'add_file=', add_file, ',',
			'out=', out_folders, ',',
			'maf_locus=', maf_locus, ',',
			'maf_finemap=', maf_finemap, ',',
			'p_hwe_finemap=', p_hwe_finemap,
			' /home/unix/dpalmer/Repositories/ldscgxe_dominance_paper/biomarker_finemapping_rerun/finemapping_template_hail_input.sh')
		)

		# Note that each of these variables: pheno, full_phenotype_file, incl_file, dom_file, add_file, and out are passed to 
		# finemapping_template_hail_input.sh and can be called using ${pheno} etc. 

		# finemapping_template_hail_input.sh takes these inputs and runs:
		# bash finemapping_hail_input.sh ${pheno} ${full_phenotype_file} ${incl_file} ${dom_file} ${add_file} ${out}

		# This generates the merged summary statistics files, together with the job scripts required to use dsub to run calculation
		# of LD and subsequent fine-mapping of the region.
	}
}

