library(data.table)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(grid)

# List the collection of biomarkers that we've run
biomarkers <- system("gsutil ls gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs-no-HWE-assumption/*", intern=TRUE)
biomarkers <- grep("both_sexes", biomarkers, value=TRUE)
biomarkers <- grep("irnt", biomarkers, value=TRUE)
biomarkers <- gsub(".*/([0-9]+).*", "\\1", biomarkers)

maf <- fread(cmd = "gsutil cat gs://ukbb-dominance/variant_AF.tsv.gz | zcat")
maf <- maf %>% mutate(variant = paste(locus, alleles ,sep=":"))
maf <- maf %>% mutate(variant = gsub("\\[", "", maf$variant))
maf <- maf %>% mutate(variant = gsub("\\]", "", maf$variant))
maf <- maf %>% mutate(variant = gsub('\\"', "", maf$variant))
maf <- maf %>% mutate(variant = gsub(',', ":", maf$variant))
maf <- maf %>% mutate(AF=gsub(".*,", "", maf$AF))
maf <- maf %>% mutate(AF=as.numeric(gsub("\\]", "", maf$AF)))
maf <- maf %>% select(variant, AF, pHWE)
maf <- maf %>% filter(AF > 0.01 & AF < 0.99) %>% filter(pHWE > 1e-6)

maf <- data.table(maf)
print(names(maf))
setkey(maf, "variant")
print(biomarkers)

biomarker_code_mapping = data.table(
	biomarker = c('Albumin', 'Alkaline phosphatase', 'Alanine aminotransferase',
	    'Apoliprotein A','Apoliprotein B', 'Aspartate aminotransferase', 'Direct bilirubin',
	    'Urea', 'Calcium', 'Cholesterol', 'Creatinine', 'C reactive protein', 'Cystatin C',
	    'Gamma glutamyltransferase', 'Glucose', 'Glycated haemoglobin', 'HDL cholesterol',
	    'IGF 1', 'LDL', 'Lipoprotein A', 'Oestradiol', 'Phosphate', 'Rheumatoid factor',
	    'SHBG', 'Total bilirubin', 'Testosterone', 'Total protein', 'Triglycerides', 'Urate',
	    'Vitamin D', 'estimated sample dilution factor'),
	biomarker_code = c('30600', '30610', '30620', '30630', '30640', '30650', '30660', '30670',
		'30680', '30690', '30700', '30710', '30720', '30730', '30740', '30750', '30760',
		'30770', '30780', '30790', '30800', '30810', '30820', '30830', '30840', '30850',
		'30860', '30870', '30880', '30890', '30897')
	)

# vp <- viewport(width = 0.35, height = 0.35, x = 0.21, y = 0.55, just=c("left", "bottom"))
# pdf(file="testing_robustness_of_HWE_GP.pdf", width=2.5, height=3)
# for (biomarker in unique(biomarkers))
# {
# 	print(biomarker)
# 	dt_HWE <- fread(cmd = paste0("gsutil cat gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs/",
# 		biomarker, "_irnt.gwas.imputed_v3.both_sexes.rerun.tsv.bgz | zcat"),
# 		select=c("variant", "dominance_pval", "dominance_tstat"), col.names=c("variant", "pval_HWE", "tstat_HWE"), key="variant")
# 	dt_HWE <- dt_HWE %>% mutate(log_pval_HWE = -((log(2) + pnorm(abs(dt_HWE[["tstat_HWE"]]), log.p=TRUE, lower.tail=FALSE)) / log(10)))

# 	dt_no_HWE <- fread(cmd = paste0("gsutil cat gs://ukb-mega-gwas-results/round2/",
# 		"dominance-biomarkers-tsvs-no-HWE-assumption/", biomarker,
# 		"_irnt.gwas.imputed_v3.both_sexes.no_HWE_assumption_GP.tsv.bgz | zcat"),
# 		select=c("variant", "dominance_pval", "dominance_tstat"), col.names=c("variant", "pval_noHWE", "tstat_noHWE"), key="variant")
# 	dt_no_HWE <- dt_no_HWE %>% mutate(log_pval_noHWE = -((log(2) + pnorm(abs(dt_no_HWE[["tstat_noHWE"]]), log.p=TRUE, lower.tail=FALSE)) / log(10)))

# 	dt <- merge(dt_HWE, dt_no_HWE)
# 	dt <- merge(maf, dt)

# 	# p <- ggplot(dt, aes(log_pval_HWE, log_pval_noHWE)) + 
# 	# geom_bin2d(bins=100) + theme_classic() + scale_fill_gradient(trans="log10", limits=c(1,1e6)) + geom_abline(intercept = 0, slope = 1, color="red")
# 	x_label <- TeX("$-\\log_{10}(\\mathit{p}_{HWE})$") 
# 	y_label <- TeX("$-\\log_{10}(\\mathit{p}_{noHWE})$")
# 	title <- biomarker_code_mapping$biomarker[which(biomarker_code_mapping$biomarker_code == biomarker)]
# 	# p <- p + labs(x=x_label, y=y_label, title=title, fill="") + theme(legend.position="none")
# 	# print(p)

# 	p <- ggplot(dt %>% filter(AF > 0.05 & AF < 0.95) %>% filter(pHWE > 1e-6), aes(log_pval_HWE, log_pval_noHWE)) + 
# 	geom_bin2d(bins=100) + theme_classic() + scale_fill_gradient(trans="log10", limits=c(1,1e6))
# 	x_label <- TeX("$-\\log_{10}(\\mathit{p}_{HWE})$") 
# 	y_label <- TeX("$-\\log_{10}(\\mathit{p}_{noHWE})$")
# 	title <- biomarker_code_mapping$biomarker[which(biomarker_code_mapping$biomarker_code == biomarker)]
# 	p <- p + labs(x=x_label, y=y_label, title=title, fill="") + theme(legend.position="none")
# 	p_and_line <- p + geom_abline(intercept = 0, slope = 1, color="red")
# 	print(p_and_line)
# 	scales <- layer_scales(p)
# 	if (min(scales$y$range$range[2], scales$x$range$range[2]) > 20) {
# 		p_inset <- p + xlim(0,10) + ylim(0,10) + labs(x="", y="", title="") +
# 		geom_abline(intercept = 0, slope = 1, color="red", size=0.2)
# 		p_inset <- p_inset +
# 		theme(
# 			panel.border = element_rect(color = "black", fill = NA, size = 0.4),
# 			panel.background = element_blank(),
# 			plot.margin = margin(t = 0,  # Top margin
# 	                             r = 0,  # Right margin
# 	                             b = 0,  # Bottom margin
# 	                             l = 0),
# 			axis.title.x = element_blank(),
# 	        axis.title.y = element_blank(),
# 	        title=element_blank(),
# 	        axis.text = element_text(size=6))
# 		print(p_inset, vp = vp)
# 	}
# }
# dev.off()

vp <- viewport(width = 0.35, height = 0.35, x = 0.21, y = 0.55, just=c("left", "bottom"))
pdf(file="testing_robustness_of_HWE_GP_0_01.pdf", width=2.5, height=3)
for (biomarker in unique(biomarkers))
{
	print(biomarker)
	dt_HWE <- fread(cmd = paste0("gsutil cat gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs/",
		biomarker, "_irnt.gwas.imputed_v3.both_sexes.rerun.tsv.bgz | zcat"),
		select=c("variant", "dominance_pval", "dominance_tstat"), col.names=c("variant", "pval_HWE", "tstat_HWE"), key="variant")
	dt_HWE <- dt_HWE %>% mutate(log_pval_HWE = -((log(2) + pnorm(abs(dt_HWE[["tstat_HWE"]]), log.p=TRUE, lower.tail=FALSE)) / log(10)))

	dt_no_HWE <- fread(cmd = paste0("gsutil cat gs://ukb-mega-gwas-results/round2/",
		"dominance-biomarkers-tsvs-no-HWE-assumption/", biomarker,
		"_irnt.gwas.imputed_v3.both_sexes.no_HWE_assumption_GP.tsv.bgz | zcat"),
		select=c("variant", "dominance_pval", "dominance_tstat"), col.names=c("variant", "pval_noHWE", "tstat_noHWE"), key="variant")
	dt_no_HWE <- dt_no_HWE %>% mutate(log_pval_noHWE = -((log(2) + pnorm(abs(dt_no_HWE[["tstat_noHWE"]]), log.p=TRUE, lower.tail=FALSE)) / log(10)))

	dt <- merge(dt_HWE, dt_no_HWE)
	dt <- merge(maf, dt)

	# p <- ggplot(dt, aes(log_pval_HWE, log_pval_noHWE)) + 
	# geom_bin2d(bins=100) + theme_classic() + scale_fill_gradient(trans="log10", limits=c(1,1e6)) + geom_abline(intercept = 0, slope = 1, color="red")
	x_label <- TeX("$-\\log_{10}(\\mathit{p}_{HWE})$") 
	y_label <- TeX("$-\\log_{10}(\\mathit{p}_{noHWE})$")
	title <- biomarker_code_mapping$biomarker[which(biomarker_code_mapping$biomarker_code == biomarker)]
	# p <- p + labs(x=x_label, y=y_label, title=title, fill="") + theme(legend.position="none")
	# print(p)

	p <- ggplot(dt %>% filter(AF > 0.01 & AF < 0.99) %>% filter(pHWE > 1e-6), aes(log_pval_HWE, log_pval_noHWE)) + 
	geom_bin2d(bins=100) + theme_classic() + scale_fill_gradient(trans="log10", limits=c(1,1e6))
	x_label <- TeX("$-\\log_{10}(\\mathit{p}_{HWE})$") 
	y_label <- TeX("$-\\log_{10}(\\mathit{p}_{noHWE})$")
	title <- biomarker_code_mapping$biomarker[which(biomarker_code_mapping$biomarker_code == biomarker)]
	p <- p + labs(x=x_label, y=y_label, title=title, fill="") + theme(legend.position="none")
	p_and_line <- p + geom_abline(intercept = 0, slope = 1, color="red")
	print(p_and_line)
	scales <- layer_scales(p)
	if (min(scales$y$range$range[2], scales$x$range$range[2]) > 20) {
		p_inset <- p + xlim(0,10) + ylim(0,10) + labs(x="", y="", title="") +
		geom_abline(intercept = 0, slope = 1, color="red", size=0.2)
		p_inset <- p_inset +
		theme(
			panel.border = element_rect(color = "black", fill = NA, size = 0.4),
			panel.background = element_blank(),
			plot.margin = margin(t = 0,  # Top margin
	                             r = 0,  # Right margin
	                             b = 0,  # Bottom margin
	                             l = 0),
			axis.title.x = element_blank(),
	        axis.title.y = element_blank(),
	        title=element_blank(),
	        axis.text = element_text(size=6))
		print(p_inset, vp = vp)
	}
}
dev.off()

