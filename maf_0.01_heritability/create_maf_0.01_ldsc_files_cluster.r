#!/usr/bin/env Rscript

library(data.table)
library(dplyr)

args = commandArgs(trailingOnly=TRUE)
batch <- as.integer(args[1])
batch_size = 20

# # Read in the set of passing phenotypes
dt_phenos <- fread("~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
phenos <- unlist(dt_phenos %>% filter(dominance_qc_no_ordinal) %>% select(phenotype))
phenos <- phenos[((batch-1)*batch_size + 1):(batch*batch_size)]
outfolder <- "/psych/genetics_data/dpalmer/UKbb/ldsc-dominance-export/"

# Extract the ldsc file

# First read in the set of SNPs to restrict to
dt_sites <- fread(
	cmd = "gsutil cat gs://ukbb-dominance/ukb_maf_auto_biallelic_maf_0.01_nomhc_snplist/b37.auto_bi_af_0.01.ukbb_gwas_qcpos.no_mhc.tsv.bgz | zcat",
	select = c("variant", "A1", "A2", "rsid"),
	col.names=c("varid", "A1", "A2", "SNP"),
	key='varid'
	)

# Additive sumstats
add_stem <- "/stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/"

# Dominance sumstats
dom_stem <- "gs://ukb-mega-gwas-results/round2/dominance-tsvs/"
dom_tsv_cloud <- system(
		"gsutil ls gs://ukb-mega-gwas-results/round2/dominance-tsvs/*both_sexes*tsv.bgz",
		intern=TRUE
		)
dom_biomarker_stem <- "gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs/"
dom_biomarker_tsv_cloud <- system("gsutil ls gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs/*rerun*tsv.bgz", intern=TRUE)

for (pheno in phenos) {

	# Then read in the additive sumstats file to grab the N
	# /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/{both_sexes}/{biomarkers}
	not_biomarker <- file.exists(paste0(add_stem, pheno, ".gwas.imputed_v3.both_sexes.tsv.bgz"))
	biomarker <- file.exists(paste0(add_stem, "biomarkers/", pheno, ".gwas.imputed_v3.both_sexes.tsv.bgz"))
	if (not_biomarker) {
		gwas_add_filepath <- paste0(add_stem, pheno, ".gwas.imputed_v3.both_sexes.tsv.bgz")
	} else if (biomarker) {
		gwas_add_filepath <- paste0(add_stem, "biomarkers/", pheno, ".gwas.imputed_v3.both_sexes.tsv.bgz")
	} else {
		print("Error: additive sumstat file not found")
		break
	}

	not_biomarker <- paste0(dom_stem, pheno, ".dominance.gwas.imputed_v3.both_sexes.tsv.bgz") %in% dom_tsv_cloud
	biomarker <- paste0(dom_biomarker_stem, pheno, ".gwas.imputed_v3.both_sexes.rerun.tsv.bgz") %in% dom_biomarker_tsv_cloud

	if (not_biomarker) {
		gwas_dom_filepath <- paste0(dom_stem, pheno, ".dominance.gwas.imputed_v3.both_sexes.tsv.bgz")
	} else if (biomarker) {
		gwas_dom_filepath <- paste0(dom_biomarker_stem, pheno, ".gwas.imputed_v3.both_sexes.rerun.tsv.bgz")
	} else {
		print("Error: dominance sumstat file not found")
		break
	}

	if (not_biomarker) {
		print("aha, not a biomarker!")
		dt_gwas_add <- fread(
			cmd = paste0("zcat ", gwas_add_filepath),
			select = c('variant', 'n_complete_samples', 'tstat'),
			col.names = c('varid', 'N', 'Z_A'),
			key = 'varid'
			)
		dt_gwas_dom <- fread(
			cmd = paste0("gsutil cat ", gwas_dom_filepath, " | zcat"),
			select = c('variant', 'dominance_tstat'),
			col.names = c('varid', 'Z_D'),
			key = 'varid'
			)
		dt_out <- merge(dt_sites, merge(dt_gwas_add, dt_gwas_dom))
		dt_out <- dt_out %>% select(SNP, A1, A2, Z_A, Z_D, N)
		fwrite(dt_out, paste0(outfolder, pheno, ".imputed_v3.ldsc_add_dom_no_hm3_filter.both_sexes.tsv.gz"), sep='\t')
	} else {
		print("A biomarker!")
	}
}

# Determine which phenotypes are missing, if any.
# List the files in the outfolder
ldsc_phenotypes <- gsub("\\.imputed.*", "", dir(outfolder))
dt_phenos <- fread("~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
phenos <- unlist(dt_phenos %>% filter(dominance_qc_no_ordinal) %>% select(phenotype))

setdiff(phenos, ldsc_phenotypes)