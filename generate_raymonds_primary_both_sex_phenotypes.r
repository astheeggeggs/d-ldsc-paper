library(data.table)
library(dplyr)

dt_original <- fread("gsutil cat gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes.tsv.bgz | gzcat")
dt_m <- fread("gsutil cat gs://ukbb-dominance/additive_dominance_gwas_results.male_including_biomarkers_corrected_nov2020.phenotypes.tsv.bgz | gzcat")
dt_f <- fread("gsutil cat gs://ukbb-dominance/additive_dominance_gwas_results.female_including_biomarkers_corrected_nov2020.phenotypes.tsv.bgz | gzcat")

# Filter out the raw continuous phenotypes and the redundant Finngen phenotypes.
dt <- dt_original %>% filter(variable_type != "continuous_raw")
dt_m <- dt_m %>% filter(variable_type != "continuous_raw")
dt_f <- dt_f %>% filter(variable_type != "continuous_raw")

redundant_finngen <- fread("raymond_finngen_redundant_phenotypes.tsv")$Phenotype_2

dt <- dt %>% filter(!(phenotype %in% redundant_finngen))
dt_m <- dt_m %>% filter(!(phenotype %in% redundant_finngen))
dt_f <- dt_f %>% filter(!(phenotype %in% redundant_finngen))

# Following these filters, in the both sexes file we have the following breakdown
			
# finngen 	501
# icd10 	633
# PHESANT 	2,922 (including 31 biomarkers, so 2,891 wihtout, matching Raymond's count)

# Determine the counts for the tests that Raymond performs

# both_sexes GWAS for phenotypes with both male and female GWAS (2714 phenotypes)
have_both <- merge(merge(dt, dt_m, by="phenotype"), dt_f, by="phenotype")
dim(have_both)
# 2,713 (which matches since we don't have the single covariate phenotype).

# both_sexes GWAS for phenotypes with neither male nor female GWAS, i.e. due to low sample size or case counts (484 phenotypes)
no_single_sex <- setdiff(dt$phenotype, union(dt_m$phenotype, dt_f$phenotype))
length(no_single_sex)
# 483 (which matches, again since we don't have the two covariate phenotypes in both_sexes, one of which isn't present in either sex).

# single-sex GWAS for phenotypes with no both_sexes GWAS (120 phenotypes)
length(setdiff(union(dt_m$phenotype, dt_f$phenotype), dt$phenotype))
# 120 (matches).

# single-sex GWAS available for phenotypes where >97% of the total sample size comes from a single sex (51 phenotypes)
table(merge(dt, dt_m, by="phenotype") %>% filter((n_non_missing.y / n_non_missing.x) > 0.97))
# binary continuous_irnt         ordinal 
#     10             262               3

male_specific_1 <- (merge(dt, dt_m, by="phenotype") %>% filter((variable_type.x != "continuous_irnt") & (n_non_missing.y / n_non_missing.x) > 0.97))$phenotype

(merge(dt, dt_m, by="phenotype") %>% filter((n_non_missing.y / n_non_missing.x) > 0.97))$phenotype
# binary continuous_irnt         ordinal 
#     30             259               8 

female_specific_1 <- (merge(dt, dt_f, by="phenotype") %>% filter((variable_type.x != "continuous_irnt") & (n_non_missing.y / n_non_missing.x) > 0.97))$phenotype

# Continuous phenotypes are not updated with the accurate numbers in the number of missing and non-missing data points (it's taken directly from the logs), so we have 
# 10 + 30 + 3 + 8 = 51.

dt <- dt %>% filter(!(phenotype %in% c(male_specific_1, female_specific_1))) 

# single-sex GWAS available for phenotypes where >99.7% of the total cases or controls comes from a single sex 
# (which is sufficient the distinguish sex-specific vs. strongly sex-biased biomedical phenotypes), excluding job codes (96 phenotypes)
male_specific_2 <- (merge(dt, dt_m, by="phenotype") %>% filter(((n_cases.y / n_cases.x) > 0.997) | ((n_controls.y / n_controls.x) > 0.997)) %>%
	filter(!grepl("Job ", description.x)))$phenotype
female_specific_2 <- (merge(dt, dt_f, by="phenotype") %>% filter(((n_cases.y / n_cases.x) > 0.997) | ((n_controls.y / n_controls.x) > 0.997)) %>% 
	filter(!grepl("Job ", description.x)))$phenotype
length(c(male_specific_2, female_specific_2))
# 96 - matches. Great!

# both_sexes GWAS for all remaining phenotypes with only one single-sex GWAS, i.e. job codes and phenotypes where
# <97% of the total samples and <99.7% of the total cases and controls are from a single sex (713 phenotypes)

# So, let's remove the male and female specific GWAS, and the collection of phenotypes that we have determined are both_sex, and 
# check that there are 713 remaining.
remainder <- setdiff(dt$phenotype, c(have_both$phenotype, no_single_sex, male_specific_1, male_specific_2, female_specific_1, female_specific_2))
# 713 phenotypes.

# Yes, everything matches. Ok, now just have to define the collection of phenotypes to restrict the both_sexes analysis to.

both_sexes_phenotypes <- c(have_both$phenotype, no_single_sex, remainder)
# 3909 - yes, this matches too.

# Now, add a true/false variable to the original file, and add in our downstream filters,and include that as a boolean too
# Easy to then merge into hail and filter with.

dt_original <- dt_original %>% 
	mutate(raymond_both_sexes_qc = (phenotype %in% both_sexes_phenotypes)) %>% 
	mutate(dominance_qc = ((raymond_both_sexes_qc) & ((!is.na(n_cases) & (n_cases > 3000) & (n_controls > 3000)) | 
                    (is.na(n_cases) & (n_non_missing > 50000))))) %>% 
	mutate(dominance_qc_no_ordinal = (dominance_qc & variable_type != "ordinal"))

# Following these additional filters, we have 1,288 phenotypes that make it into the heritability analysis if we include ordinal variables,
# and 1,060 if we don't.
fwrite(
	dt_original,
	file="additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv",
	sep="\t"
	)
