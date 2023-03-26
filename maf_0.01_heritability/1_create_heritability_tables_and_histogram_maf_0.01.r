library(data.table)
library(dplyr)
library(ggplot2)

# Additive and dominance
# Weights and Ref the same - summed over only HM3 SNPs
h2 <- fread("h2_estimates/add_dom_maf_0.01.h2")

# Remove all of the case control phenotypes that have less than 3000 cases.
phenotype_filters <- fread("../../ldscgxe_dominance_paper/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv")
h2 <- merge(h2, phenotype_filters, by=c("phenotype"))

# Restrict to phenotypes with sufficient samples, and that pass Raymonds checks. 
# Need to make sure that this collection of phenotypes is correct - need to include analysis of biomarkers.

get_liab <- function(K, h2_obs) {
	return(ifelse(is.na(K), h2_obs, h2_obs * K * (1 - K) / (dnorm(qnorm(K))^2)))
}

h2$K <- h2$n_cases/(h2$n_cases + h2$n_controls)

for(i in 1:nrow(h2)) {

	# Liability scale heritability estimates.
	h2$h2_liability_A[i] <- get_liab(h2$K[i], h2$h2_observed_A[i])
	h2$h2_liability_D[i] <- get_liab(h2$K[i], h2$h2_observed_D[i])

	# Liability scale standard errors.
	h2$h2_liability_se_A[i] <- get_liab(h2$K[i], h2$h2_observed_se_A[i])
	h2$h2_liability_se_D[i] <- get_liab(h2$K[i], h2$h2_observed_se_D[i])

}

# Write the full h2 data to disk - this is then used to generate the mean estimate and 
# 95% confidence intervals that we want for the abstract.
fwrite(h2, file="~/Repositories/ldscgxe_dominance_paper/h2_table_maf_0.01.tsv", sep="\t")
h2 <- fread("~/Repositories/ldscgxe_dominance_paper/h2_table_maf_0.01.tsv")

library(IsoplotR)
york(cbind(h2$h2_liability_A,h2$h2_liability_se_A,h2$h2_liability_D,h2$h2_liability_se_D, cor(h2$h2_liability_se_A, h2$h2_liability_se_D)))

# library(yorkregression)
# york_fit <- york(x=h2$h2_liability_A, y=h2$h2_liability_D, sd_x=h2$h2_liability_se_A, sd_y=h2$h2_liability_se_D, r_xy_errors=cor(h2$h2_liability_se_A, h2$h2_liability_se_D))

lm(h2_liability_D ~ h2_liability_A + 0, data=h2)
lm(h2_liability_D ~ h2_liability_A + 0, data=h2 %>% filter(variable_type == "continuous_irnt"))

pdf(file='~/Repositories/ldscgxe_dominance_paper/maf_0.01_heritability/heritability_results_summaries/h2_histograms.pdf', width=4, height=4)
p <- ggplot(h2) + 
	geom_histogram(aes(x = h2_liability_A),
		alpha=0.6, fill ="cornflowerblue", binwidth=0.01, position="dodge") + 
	geom_histogram(aes(x = h2_liability_D),
		alpha=0.6, fill ="indianred3", binwidth=0.01, position="dodge") + 
	xlab("Heritability estimate") + ylab("Frequency") + theme_classic() +
	geom_vline(mapping=aes(xintercept=mean(h2_liability_A)), color='grey') +
	geom_vline(mapping=aes(xintercept=mean(h2_liability_D)), color='grey')
print(p)
dev.off()

pdf(file='~/Repositories/ldscgxe_dominance_paper/maf_0.01_heritability/heritability_results_summaries/h2_scatter.pdf', width=4, height=4)
p <- ggplot(h2, aes(x = h2_liability_A, y = h2_liability_D)) + 
	geom_errorbarh(aes(xmin = h2_liability_A - 1.96 * h2_liability_se_A, xmax = h2_liability_A + 1.96 * h2_liability_se_A), color="cornflowerblue") +
	geom_errorbar(aes(ymin = h2_liability_D - 1.96 * h2_liability_se_D, ymax = h2_liability_D + 1.96 * h2_liability_se_D), color="indianred3") + 
	geom_point(color="grey", cex=0.5) + 
	xlab("Additive heritability estimate") + ylab("Dominance heritability estimate") + theme_classic() + 
	geom_abline(slope=york_fit$coefficient[2,1], intercept=york_fit$coefficient[1,1], lwd=1) + 
	coord_cartesian(xlim=c(0,1), ylim = c(-0.2, NA))
print(p)
dev.off()


# Get the table of phenotypes at the top end for both additive and dominance results.

# Automate to create a latex table.
# Create table from the h2 output

# Include the p-value, standard error and h2 estimate.

# Filter down to the curated phenotype list.
h2 <- h2 %>% filter(dominance_qc_no_ordinal)

h2_D <- h2 %>% select(description, h2_liability_D, h2_liability_se_D,  h2_p_D) %>% arrange(h2_p_D) %>% head(100) %>% 
	mutate(
	h2_liability_D = signif(h2_liability_D, 3),
	h2_liability_se_D = signif(h2_liability_se_D, 3),
	h2_p_D = signif(h2_p_D, 3)
	)
	
h2_A <- h2 %>% select(description, h2_liability_A, h2_liability_se_A,  h2_p_A) %>% arrange(h2_p_A) %>% head(100) %>% 
	mutate(
	h2_liability_A = signif(h2_liability_A, 3),
	h2_liability_se_A = signif(h2_liability_se_A, 3),
	h2_p_A = signif(h2_p_A, 3)
	)

h2_D <- h2_D %>% rename("Phenotype"= description,
	"$h^2_D$"= h2_liability_D,
	"SE $h^2_D$"= h2_liability_se_D,
	"$p$ $h^2_D$"= h2_p_D)

h2_A <- h2_A %>% rename("Phenotype"= description,
	"$h^2_A$"= h2_liability_A,
	"SE $h^2_A$"= h2_liability_se_A,
	"$p$ $h^2_A$"= h2_p_A)

con <- file("h2_results_A_maf_0.01.tex", "w")

writeLines(c("\\begin{table}[]", 
    "\\centering",
    "\\begin{tabular}{l l l l l l l l l}",
    paste(paste(names(h2_A), collapse=" & "), "\\\\\\hline")), con = con)

for(i in 1:nrow(h2_A))
{
	writeLines(paste(paste(h2_A[i,], collapse=" & "), "\\\\"), con = con)
}
writeLines(c("\\hline", "\\end{tabular}", "\\caption{}", "\\end{table}"), con=con)

close(con)


con <- file("h2_results_D_maf_0.01.tex", "w")

writeLines(c("\\begin{table}[]", 
    "\\centering",
    "\\begin{tabular}{l l l l l l l l l}",
    paste(paste(names(h2_D), collapse=" & "), "\\\\\\hline")), con = con)

for(i in 1:nrow(h2_D))
{
	writeLines(paste(paste(h2_D[i,], collapse=" & "), "\\\\"), con = con)
}

writeLines(c("\\hline", "\\end{tabular}", "\\caption{}", "\\end{table}"), con=con)

close(con)
