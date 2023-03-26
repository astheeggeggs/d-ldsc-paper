library(data.table)
library(dplyr)
library(plotly)
library(processx)
library(ggplot2)

# DEV: Sanity checked that beta_A is on the standardised scale - This is done at the end of the get_top_dominance_and_additive_loci_for_all_phenotypes.py, 
# I also double checked it in hail.

read_and_restrict <- function(gs_location,
  dt_pheno_location="~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv",
  add_hit_restrict="pval_A", dom_hit_restrict="pval_D", type_name="both", sig=5e-8, maf_0.01=FALSE) {
  dt <- fread(cmd=paste0('gsutil cat ', gs_location, ' | gzcat'))

  # Merge in the phenotype data that we have and restrict to the filtered phenotypes.
  dt_pheno <- fread(dt_pheno_location)
  dt_pheno <- dt_pheno %>% filter(raymond_both_sexes_qc) %>% select(phenotype, dominance_qc, dominance_qc_no_ordinal)

  dt <- merge(dt, dt_pheno)
  if (maf_0.01) {
    dt <- dt %>% filter((n_cases > 75000 & n_controls > 75000) | (is.na(n_cases) & n_non_missing > 75000))
  }

  dt$variable_type[which(dt$variable_type == "binary")] <- "Binary"
  dt$variable_type[which(dt$variable_type == "categorical")] <- "Binary"
  dt$variable_type[which(dt$variable_type == "continuous_irnt")] <- "Continuous (IRNT)"
  dt$variable_type[which(dt$variable_type == "continuous_raw")] <- "Continuous (raw)"
  dt$variable_type[which(dt$variable_type == "ordinal")] <- "Ordinal"
  dt$variable_type <- as.factor(dt$variable_type)

  pheno_description <- paste0(dt$phenotype, ' : ', dt$description)
  dt_filter <- dt %>% filter(dominance_qc_no_ordinal)

  V_A <- rowSums(select(dt_filter, starts_with("beta_A")) ** 2)
  V_D <- rowSums(select(dt_filter, starts_with("beta_D")) ** 2)

  dt_V <- data.table(type=type_name, V_A=V_A, V_D=V_D, just_hits=FALSE) %>% filter(!((V_A == 0) & (V_D == 0)))

  # Want to create a plot of just the hits...
  V_A_hits <- rowSums((select(dt_filter, starts_with("beta_A")) ** 2) * ((select(dt_filter, starts_with(add_hit_restrict))) < sig))
  V_D_hits <- rowSums((select(dt_filter, starts_with("beta_D")) ** 2) * ((select(dt_filter, starts_with(dom_hit_restrict))) < sig))

  dt_V <- merge(dt_V, data.table(type=type_name, V_A=V_A_hits, V_D=V_D_hits, just_hits=TRUE) %>% filter(!((V_A == 0) | (V_D == 0))), all=TRUE)

  return(dt_V)
}

# dt_V <- read_and_restrict(gs_location="gs://ukbb-dominance/additive.dominance.gwas.imputed_v3.both_sexes_including_biomarkers_corrected_jan2020_top_hits_no_X.tsv.bgz", sig=1e-6)
dt_V <- read_and_restrict(gs_location="gs://ukbb-dominance/additive.dominance.gwas.imputed_v3.both_sexes_including_nov2022_biomarkers_corrected_jan2020_top_hits_no_X_maf_0.05.tsv.bgz", sig=1e-6)
median_lines <- dt_V %>% group_by(just_hits) %>% summarise(z = median(V_A/V_D))

pdf("~/Repositories/ldscgxe_dominance_paper/figure_3/hist_of_V_A_V_D_facet_nov2022_biomarkers_maf_0.05.pdf", width=7, height=3.5)
p <- ggplot(dt_V, aes(x=V_A/V_D)) + geom_histogram(fill="cornflowerblue", color="grey") + scale_x_log10() + 
  ylab("Frequency") + theme_classic() + facet_wrap(vars(just_hits))
p <- p + geom_vline(aes(xintercept=z), data=median_lines, size=1.5, color="red")
print(p)
dev.off()

dt_V <- read_and_restrict(
  gs_location="gs://ukbb-dominance/additive.dominance.gwas.imputed_v3.both_sexes_including_nov2022_biomarkers_corrected_jan2020_top_hits_no_X_maf_0.01.tsv.bgz",
  sig=1e-6, maf_0.01=TRUE)
median_lines <- dt_V %>% group_by(just_hits) %>% summarise(z = median(V_A/V_D))

pdf("~/Repositories/ldscgxe_dominance_paper/figure_3/hist_of_V_A_V_D_facet_nov2022_biomarkers_maf_0.01.pdf", width=7, height=3.5)
p <- ggplot(dt_V, aes(x=V_A/V_D)) + geom_histogram(fill="cornflowerblue", color="grey") + scale_x_log10() + 
  ylab("Frequency") + theme_classic() + facet_wrap(vars(just_hits))
p <- p + geom_vline(aes(xintercept=z), data=median_lines, size=1.5, color="red")
print(p)
dev.off()