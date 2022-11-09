library(data.table)
library(dplyr)

# Run it
path_to_output <- "/stanley/genetics_storage/analysis/ukbb_dominance/sumstats-dominance-export/both_sexes/"

# First get the filelist
gs_files <- system("gsutil ls gs://ukb-mega-gwas-results/round2/dominance-tsvs/*both_sexes*tsv.bgz", intern=TRUE)
dt_gs <- data.table(gs_folder=gs_files)
dt_gs <- dt_gs %>% mutate(phenotype = gsub(".*dominance-tsvs/(.*).dominance.*", "\\1", dt_gs$gs_folder))
dt_gs <- data.table(dt_gs)
setkey(dt_gs, "phenotype")
# Then, filter to the set of files that we want to download from the cloud
dt_phenos <- fread("~/Repositories/ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
phenos <- dt_phenos %>% filter(dominance_qc_no_ordinal) %>% select(phenotype)
phenos <- data.table(phenos)
setkey(phenos, "phenotype")

dt_gs_filtered <- merge(phenos, dt_gs)

# Then, loop over each of them, submitting the job to the cluster

for (gs_file in dt_gs_filtered$gs_folder[900:length(dt_gs_filtered$gs_folder)]) {
	system(paste0("qsub -v gs_file=", gs_file, " move_dom_sumstats_from_cluster_template.sh"))
}
