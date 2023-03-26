library(dplyr)
library(data.table)
source("~/Repositories/BRaVa_curation/QC_UKBB_200k_BMRC/utils/pretty_plotting.r")
library(ggridges)
library(latex2exp)

h2 <- fread("~/Repositories/ldscgxe_dominance_paper/h2_table.tsv") %>% filter(dominance_qc_no_ordinal)

h2_meta <- h2 %>% select(
	phenotype, n_missing, n_non_missing, PHESANT_transformation,
	n_cases, n_controls, variable_type, source, notes
	)
setkey(h2_meta, "phenotype")
phenotypes <- unique(h2$phenotype)

# Read in all the estimates
dt <- data.table()
for (a in c("-1.0", "-0.5", "-0.25", "0.0", "0.25", "0.5", "1.0"))
{
	dt_add_tmp <- fread(paste0("h2_estimates/add.alpha.", a, ".h2")) %>% mutate(alpha=a, type="add")
	dt_dom_tmp <- fread(paste0("h2_estimates/dom.alpha.", a, ".h2")) %>% mutate(alpha=a, type="dom")
	dt_add_tmp <- data.table(dt_add_tmp)
	dt_dom_tmp <- data.table(dt_dom_tmp)
	dt_tmp <- rbind(dt_add_tmp, dt_dom_tmp) %>% filter(phenotype %in% phenotypes)
	dt <- rbind(dt, dt_tmp)
}

dt <- dt %>% mutate(
	h2_p = -log10(h2_p),
	alpha = as.factor(as.numeric(dt$alpha)-1)
	)

ribbon_p <- 0.95
dt_plot <- data.table(
	dt %>% arrange(type, alpha, desc(h2_p)) %>% select(-c(type,alpha)),
	dt %>% group_by(type, alpha) %>% arrange(desc(h2_p)) %>% 
		summarize(
			Pvalue_expected = -log10(seq(1, n())/(n() + 1)),
			clower = -log10(qbeta(p = (1 - ribbon_p) / 2, shape2 = n():1, shape1 = 1:n())),
			cupper = -log10(qbeta(p = (1 + ribbon_p) / 2, shape2 = n():1, shape1 = 1:n()))
			)
		)


create_pretty_qq_plot(
	plot_title="Additive",
	cex_labels=2,
	dt_plot %>% filter(type=="add"), aes(x=Pvalue_expected, y=h2_p, color=alpha),
	save_figure=TRUE,
	x_label=TeX("Expected $(-\\log_{10}\\mathit{P})$"),
	y_label=TeX("Observed $(-\\log_{10}\\mathit{P})$"),
	key_cols=NULL,
	aes_ribbon = aes(ymin=clower, ymax=cupper),
	width=120, height=90,
	file="plots/additive_heritability_qq_plot",
	raster=FALSE
)

create_pretty_qq_plot(
	plot_title="Dominance",
	cex_labels=2,
	dt_plot %>% filter(type=="dom"), aes(x=Pvalue_expected, y=h2_p, color=alpha),
	save_figure=TRUE,
	x_label=TeX("Expected $(-\\log_{10}\\mathit{P})$"),
	y_label=TeX("Observed $(-\\log_{10}\\mathit{P})$"),
	key_cols=c("alpha"),
	aes_ribbon = aes(ymin=clower, ymax=cupper),
	width=120, height=90,
	file="plots/dominance_heritability_qq_plot",
	raster=FALSE
)

get_liab <- function(K, h2_obs) {
	return(ifelse(is.na(K), h2_obs, h2_obs * K * (1 - K) / (dnorm(qnorm(K))^2)))
}

setkey(dt, "phenotype")
dt <- merge(dt, h2_meta)

dt <- dt %>% mutate(
	K = n_cases/(n_cases + n_controls),
	h2_liability = get_liab(K, h2_observed),
	h2_liability_se = get_liab(K, h2_observed_se)
)

pdf(file='plots/h2_histograms.pdf', width=4, height=3)
p <- ggplot(dt %>% filter(type=="add")) + 
	geom_density_ridges(aes(x = h2_liability, y=alpha, color=alpha, fill=alpha),
		alpha=0.6) + labs(x="Additive heritability estimate") + theme_classic() +
        scale_color_d3('category20') + scale_fill_d3('category20')
print(p)
p <- ggplot(dt %>% filter(type=="dom")) + 
	geom_density_ridges(aes(x = h2_liability, y=alpha, color=alpha, fill=alpha),
		alpha=0.6) + labs(x="Dominance heritability estimate") + theme_classic() +
        scale_color_d3('category20') + scale_fill_d3('category20')
print(p)
dev.off()

# Create a correlation plot of the heritabilities across the alphas
library(corrplot)
dt_cor <- dcast(dt, phenotype ~ alpha + type, value.var="h2_liability")
cor_add <- cor(dt_cor %>% select(grep("add", names(dt_cor), value=TRUE)), use="complete.obs")
cor_dom <- cor(dt_cor %>% select(grep("dom", names(dt_cor), value=TRUE)), use="complete.obs")
rownames(cor_add) <- colnames(cor_add) <- gsub("_add", "", colnames(cor_add))
rownames(cor_dom) <- colnames(cor_dom) <- gsub("_dom", "", colnames(cor_dom))

pdf(file='plots/h2_corr.pdf', width=5, height=4)
corrplot(cor_add, type = "upper", order = "original", method="number",
         tl.col = "black")
corrplot(cor_dom, type = "upper", order = "original", method="number",
         tl.col = "black")
dev.off()

# Plot the LD-scores against each other
dt_ld_list <- list()
i <- 1
for (a in c("-1.0", "-0.5", "-0.25", "0.0", "0.25", "0.5", "1.0")) {
	print(a)
	for (chr in seq(1,22)) {
		print(chr)
		dt_ld_list[[i]] <- fread(paste0("ldscores/1000G.EUR.QC.alpha.", a, ".chr", chr, ".ldscore")) %>% mutate(alpha=a)
		names(dt_ld_list[[i]]) <- c("CHR", "SNP", "BP", "L2_A", "L2_D", "alpha")
		dt_ld_list[[i]] <- melt(dt_ld_list[[i]], id.vars=c("SNP", "alpha"), measure.vars=c("L2_A", "L2_D"))
		i <- i+1
	}
}
dt_ld_list <- rbindlist(dt_ld_list)
dt_ld_list$alpha = as.factor(as.numeric(dt_ld_list$alpha)-1)

pdf(file='plots/ld_score_dists.pdf', width=4, height=3)
p <- ggplot(dt_ld_list) + 
	geom_density_ridges(aes(x = value, y=alpha, color=variable, fill=variable),
		alpha=0.6) + labs(x="LD-score") + theme_classic() +
        scale_color_manual(values=c("cornflowerblue", "indianred3")) + scale_fill_manual(values=c("cornflowerblue", "indianred3")) + xlim(c(-10, 100))
p <- p + geom_vline(xintercept = 1, color = "grey10", size=1)
print(p)
dev.off()

