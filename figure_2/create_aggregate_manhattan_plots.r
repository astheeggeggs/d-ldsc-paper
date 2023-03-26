library(data.table)
library(dplyr)
library(ggplot2)
library(ggrastr)

source("lightweight_manhattan.r")

dt_dict <- fread("~/Repositories/PHESANT/variable-info/Data_Dictionary_Showcase_Prioritized_list_for Pfizer_and Biogen_v3.0.tsv", sep='\t', key='FieldID')

# system("gsutil cp gs://ukbb-dominance/aggregate_manhattan_no_ord_maf_0_05.tsv.bgz .")
# system("gsutil cp gs://ukbb-dominance/aggregate_manhattan_maf_0_05.tsv.bgz .")
system("gsutil cp gs://ukbb-dominance/aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz .")
system("gsutil cp gs://ukbb-dominance/aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz .")

# Comment out - significantly faster if we download first and read in.
# commands <- c(
# 	# "gsutil cat gs://ukbb-dominance/aggregate_manhattan_maf_0_05.tsv.bgz | gzcat",
# 	# "gsutil cat gs://ukbb-dominance/aggregate_manhattan_no_ord_maf_0_05.tsv.bgz | gzcat")

# commands <- c(
# 	"gzcat aggregate_manhattan_maf_0_05.tsv.bgz",
# 	"gzcat aggregate_manhattan_no_ord_maf_0_05.tsv.bgz"
# 	)

commands <- c(
	"gzcat aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz",
	"gzcat aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz"
	)

for(command in commands) {
	dt <- fread(command)

	dt <- dt %>% mutate(extreme_phenotype_dominance_to_match = gsub("\\[\"(.*)\"\\]", "\\1", extreme_phenotype_dominance),
		extreme_phenotype_additive_to_match = gsub("\\[\"(.*)\"\\]", "\\1", extreme_phenotype_additive))
	dt <- dt %>% mutate(extreme_phenotype_dominance_to_match = gsub("_[0-9]+", "", extreme_phenotype_dominance_to_match),
		extreme_phenotype_additive_to_match = gsub("_[0-9]+", "", extreme_phenotype_additive_to_match))
	dt <- dt %>% mutate(extreme_phenotype_dominance_to_match = gsub("_irnt", "", extreme_phenotype_dominance_to_match),
		extreme_phenotype_additive_to_match = gsub("_irnt", "", extreme_phenotype_additive_to_match))

	# Create a new column that's the category of the phenotype and use that for colouring the manhattan plot.
	dt_dict <- dt_dict %>% mutate(curated_category = ifelse(sapply(strsplit(dt_dict$Path, split=" > "), `[`, 2) == "Touchscreen", paste(sapply(strsplit(dt_dict$Path, split=" > "), `[`, 3), "(touchscreen)"),
		ifelse(sapply(strsplit(dt_dict$Path, split=" > "), `[`, 2) == "Verbal interview", paste(sapply(strsplit(dt_dict$Path, split=" > "), `[`, 3), "(verbal interview)"),
		sapply(strsplit(dt_dict$Path, split=" > "), `[`, 2))))
	dt_dict <- data.table(dt_dict)
	setkey(dt_dict, "FieldID")

	dt <- dt %>% mutate(
		dominance_category = ifelse(
			is.na(dt_dict[extreme_phenotype_dominance_to_match]$curated_category),
				ifelse(nchar(extreme_phenotype_dominance_to_match) == 3, 
					"ICD10 code",
					ifelse(grepl("^[A-Z]", extreme_phenotype_dominance_to_match),
						"FinnGen encoding", "Biomarkers"
					)
				),
			dt_dict[extreme_phenotype_dominance_to_match]$curated_category),
		additive_category = ifelse(
			is.na(dt_dict[extreme_phenotype_additive_to_match]$curated_category),
				ifelse(nchar(extreme_phenotype_additive_to_match) == 3, 
					"ICD10 code",
					ifelse(grepl("^[A-Z]", extreme_phenotype_additive_to_match),
						"FinnGen encoding", "Biomarkers"
					)
				),
			dt_dict[extreme_phenotype_additive_to_match]$curated_category)
		)

	to_write <- gsub(".*(aggregate.*tsv).*", "\\1", command)
	fwrite(dt, sep='\t', file=to_write)
	system(paste("gzip", to_write))
	system(paste0("mv ", to_write, ".gz ", to_write, ".bgz"))
}

# make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance", colour_col="dominance_category", transform=loglogtrans, breaks=c(seq(0, 30, by=5), 10^seq(3,4)))
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive", colour_col="additive_category", transform=loglogtrans2, breaks=c(seq(0, 300, by=50), 10^seq(3,4)))
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_transform", colour_col="dominance_category")
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_transform", colour_col="additive_category")

# make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_ord", colour_col="dominance_category", transform=loglogtrans, breaks=c(seq(0, 30, by=5), 10^seq(3,4)))
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_ord", colour_col="additive_category", transform=loglogtrans2, breaks=c(seq(0, 300, by=50), 10^seq(3,4)))
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_transform_no_ord", colour_col="dominance_category")
# make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_transform_no_ord", colour_col="additive_category")

make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance", colour_col="dominance_category", transform=loglogtrans, breaks=c(seq(0, 30, by=5), 10^seq(3,4)))
make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive", colour_col="additive_category", transform=loglogtrans2, breaks=c(seq(0, 300, by=50), 10^seq(3,4)))
make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_transform", colour_col="dominance_category")
make_ukbb_summary_manhattan_plot("aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_transform", colour_col="additive_category")

make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_ord", colour_col="dominance_category", transform=loglogtrans, breaks=c(seq(0, 30, by=5), 10^seq(3,4)))
make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_ord", colour_col="additive_category", transform=loglogtrans2, breaks=c(seq(0, 300, by=50), 10^seq(3,4)))
make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_dominance_tstat", plot_title="Dominance", plot_name="meta_dominance_no_transform_no_ord", colour_col="dominance_category")
make_ukbb_summary_manhattan_plot("aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz", tstat_col="extreme_additive_tstat", plot_title="Additive", plot_name="meta_additive_no_transform_no_ord", colour_col="additive_category")

# Create colour scheme.

# Create plots - annotating with the largest hits.
# Merge with the Annovar to get the phenotypes and the closest gene for the top hits (>1e-30)

# To label the hits.
# Grab this file from the Broad cluster - note that it excludes the HLA region as the LD extends so far.
# scp dpalmer@login01:/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping_biomarkers_rerun//lead_dominance_variant_maf_0.05_p_hwe_1e-06_maf_finemap_0.01.tsv .
dt <- fread("lead_dominance_variant_maf_0.05_p_hwe_1e-06_maf_finemap_0.01.tsv")
dt <- dt %>% filter(log_p_dominance < -30)

dt <- dt %>% mutate(variant=paste(chromosome, position, allele1, allele2, sep=":"))
setkey(dt, "phenotype")
dt_pheno <- fread("~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
setkey(dt_pheno, "phenotype")
dt <- merge(dt, dt_pheno) %>% filter(dominance_qc_no_ordinal)
dt <- dt %>% filter(!duplicated(rsid)) %>% select(phenotype, variant, rsid, log_p_additive, log_p_dominance, Gene.refGene, phenotype_description, maf)

dt <- dt %>% select(variant,  Gene.refGene)
setkey(dt, "variant")
