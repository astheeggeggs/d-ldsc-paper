library(data.table)
library(dplyr)

# Create the .incl files for each of the different significant loci files

# First, check to ensure that all of the phenotype files on the cluster are the same as those on the cloud.

# Full phenotypes files
biomarkers_cloud <- "gs://ukb31063/ukb31063.biomarkers.tsv.bgz"
phesant_cloud <- "gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz"
ICD10_cloud <- "gs://ukb31063/ukb31063.ICD10_phenotypes.both_sexes.tsv.bgz"
finngen_cloud <- "gs://ukb31063/ukb31063.FINNGEN_phenotypes.both_sexes.tsv.bgz"

full_phenotype_files_cloud <- c(phesant_cloud, ICD10_cloud, finngen_cloud, biomarkers_cloud)

# Updated to the archived location
biomarkers_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.biomarkers.tsv.bgz"
phesant_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz"
ICD10_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.ICD10_phenotypes.both_sexes.tsv.bgz"
finngen_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.FINNGEN_phenotypes.both_sexes.tsv.bgz"

full_phenotype_files_cluster <- c(phesant_cluster, ICD10_cluster, finngen_cluster, biomarkers_cluster)

# Summary of phenotype files
phesant_summary_cloud <- "gs://ukb31063/ukb31063.PHESANT_January_2019.both_sexes.summary.tsv"
ICD10_summary_cloud <- "gs://ukb31063/ukb31063.ICD10_phenotypes.both_sexes.summary.tsv"
finngen_summary_cloud <- "gs://ukb31063/ukb31063.FINNGEN_phenotypes.both_sexes.summary.tsv"

phenotype_summaries_cloud <- c(phesant_summary_cloud, ICD10_summary_cloud, finngen_summary_cloud)

phesant_summary_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.PHESANT_January_2019.both_sexes.summary.tsv"
ICD10_summary_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.ICD10_phenotypes.both_sexes.summary.tsv"
finngen_summary_cluster <- "/stanley/genetics_storage/ukb31063_archive/ukb31063/ukb31063.FINNGEN_phenotypes.both_sexes.summary.tsv"

phenotype_summaries_cluster <- c(phesant_summary_cluster, ICD10_summary_cluster, finngen_summary_cluster)

# This is a new location
out_folders <- c(
	"/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/inputs/incl_files/phesant",
	"/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/inputs/incl_files/icd10",
	"/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/inputs/incl_files/finngen",
	"/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun/inputs/incl_files/biomarkers"
	)

for (folder in out_folders) {
	system(paste0("mkdir -p ", folder))
}

awk_script <- "~/Repositories/ldscgxe/dominance_fine_mapping/incl_file_preparation.awk"

for(i in 1:length(full_phenotype_files_cloud)) {
	if(gsub(" .*", "", system(paste0("gsutil cat ", full_phenotype_files_cloud[i], " | md5sum"), intern=TRUE)[[1]]) != gsub(" .*", "", system(paste("md5sum", full_phenotype_files_cluster[i]), intern=TRUE)[[1]])) {
		cat("Error, md5checksums do not match!\n")
	} else {
		cat("Pass, md5checksums match.\n")
	}
}

for(i in 1:length(phenotype_summaries_cloud)) {
	if(gsub(" .*", "", system(paste0("gsutil cat ", phenotype_summaries_cloud[i], " | md5sum"), intern=TRUE)[[1]]) != gsub(" .*", "", system(paste("md5sum", phenotype_summaries_cluster[i]), intern=TRUE)[[1]])) {
		cat("Error, md5checksums do not match!\n")
	} else {
		cat("Pass, md5checksums match.\n")
	}
}

full_phenotype_files_cluster[4] <- "/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/ukb30163.restricted_and_phesant_formatted.biomarkers.tsv.gz"

# phenotype_dt <- fread(cmd = "gsutil cat gs://ukbb-dominance/phenotypes_with_significant_dominance_hit_no_cts_raw.tsv.bgz | zcat")
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
	# Find out which full phenotype file the current phenotype is in.
	k <- 0
	for (j in 4:length(full_phenotype_files_cluster)) { # -1 because the biomarkers behave slightly differently.
		# Check to see if this phenotype is in this full phenotype file.
		file_phenotypes <- strsplit(system(paste("zcat", full_phenotype_files_cluster[j], "| head -1"), intern=TRUE)[[1]], split='\t')[[1]]
		# Different behavior for the biomarkers.
		if(j==4) {
			file_phenotypes <- paste0(file_phenotypes, "_irnt")
		}

		if(phenotypes[i] %in% file_phenotypes) {
			out_file_tmp <- paste0(out_folders[j], "/", phenotypes[i], "_tmp.incl")
			out_file <- paste0(out_folders[j], "/", phenotypes[i], ".incl")
			cat("Current phenotype:", phenotypes[i], "\n")
			cat("Match found in:", full_phenotype_files_cluster[j], "\n")
			system(paste0("zcat ",  full_phenotype_files_cluster[j], " | awk -v col=", ifelse(j==4, gsub("_irnt", "", phenotypes[i]), phenotypes[i]), " -f ", awk_script, "> ", out_file_tmp))

			if(j!=4) { # Different behavior for the biomarkers.
				# Run a test to make sure the counts agree with the summary file.
				n_lines <- strsplit(system(paste0("wc -l ", out_file_tmp), intern=TRUE), split=" ")[[1]][1]
				check <- strsplit(system(paste0("awk 'BEGIN { FS = \"\t\" }; { print $1\"\t\"$3 }' ", phenotype_summaries_cluster[j], " | awk -v phe=", phenotypes[i], " 'BEGIN { FS = \"\t\" } ; $1==phe'"), intern=TRUE), split='\t')[[1]][2]
				
				if(n_lines != check) {
					cat("Error, problem with phenotype:", phenotypes[i], "\n")
				} else {
					cat("Pass, counts match for phenotype:", phenotypes[i], n_lines, "\n")
				}
			} else {
				# Run a test to make sure the counts agree with the summary file.
				n_lines <- strsplit(system(paste0("wc -l ", out_file_tmp), intern=TRUE), split=" ")[[1]][1]
				check <- strsplit(system(paste0("gsutil cat gs://ukbb-dominance/biomarkers.both_sexes.tsv | awk 'BEGIN { FS = \"\t\" }; { print $1\"\t\"$5 }' | awk -v phe=", phenotypes[i], " 'BEGIN { FS = \"\t\" } ; $1==phe'"), intern=TRUE), split='\t')[[1]][2]
				
				if(n_lines != check) {
					cat("Error, problem with phenotype:", phenotypes[i], "\n")
				} else {
					cat("Pass, counts match for phenotype:", phenotypes[i], n_lines, "\n")
				}
			}

			if(!file.exists(out_file)) {
				# If the file doesn't exist, remove the _tmp extension.
				cat("This phenotype:", phenotypes[i], "; incl file does not exist, removing _tmp extension.\n")
				system(paste("mv", out_file_tmp, out_file))
			} else {
				# If the file does exist, check to see if the newly created file is the same as the old one.
				if(gsub(" .*", "", system(paste("md5sum", out_file_tmp), intern=TRUE)[[1]]) == gsub(" .*", "", system(paste("md5sum", out_file), intern=TRUE)[[1]])) {
					# If it is, remove the one with the _tmp extension.
					cat("This phenotype:", phenotypes[i], "; .incl file exists and is identical to previous. Removing _tmp.incl file.\n")
					system(paste("rm", out_file_tmp))
				} else {
					cat("Error, there is a difference in the samples determined for .incl and _tmp.incl for phenotype:", phenotypes[i], "\n")
				}
			}
			k <- k+1
		}
		
	}

	if(k == 0) cat("Error, none of the phenotype files found a match for phenotype:", phenotypes[i], ".\n")
	if(k > 1) cat("Error, multiple matches for this phenotype:", phenotypes[i], "across full phenotype files.\n")

}

