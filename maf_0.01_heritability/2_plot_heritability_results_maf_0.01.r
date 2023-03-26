library(data.table)
library(dplyr)
library(ggplot2)
library(viridis)

qqconf <- function(n, conf.points=1000, conf.col="gray", conf.alpha=.05)
{
	conf.points = min(conf.points, n-1)
	mpts <- matrix(nrow=conf.points*2, ncol=2)

	seq_points <- seq(from=1, to=conf.points)
	mpts[seq_points,1] <- -log10((seq_points-0.5)/n)
	mpts[seq_points,2] <- -log10(qbeta(1-conf.alpha/2, seq_points, n-seq_points))
	mpts[conf.points*2+1-seq_points,1] <- -log10((seq_points-0.5)/n)
	mpts[conf.points*2+1-seq_points,2] <- -log10(qbeta(conf.alpha/2, seq_points, n-seq_points))

	polygon(x=mpts[,1], y=mpts[,2], border=FALSE, col='grey90')
}

h2 <- fread("~/Repositories/ldscgxe_dominance_paper/h2_table_maf_0.01.tsv")

cex_point=0.5
line_width=0.5

pdf("~/Repositories/ldscgxe_dominance_paper/maf_0.01_heritability/heritability_results_summaries/additive_dominance_maf_0.01_qq_plots.pdf", width=9, height=4.5)
par(mfrow=c(1,2))

plot(sort(-log10(seq(0,1,length.out=(length(h2$h2_p_A)+1))[-1])),
	sort(-log10(h2$h2_p_A)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")), ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Additive', col='cornflowerblue', cex=cex_point)
qqconf(length(h2$h2_p_A)+1)
points(sort(-log10(seq(0,1,length.out=(length(h2$h2_p_A)+1))[-1])),
	sort(-log10(h2$h2_p_A)),
	pch=16, col='cornflowerblue', cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)

plot(sort(-log10(seq(0,1,length.out=(length(h2$h2_p_D)+1))[-1])),
	sort(-log10(h2$h2_p_D)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")), ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Dominance', col='cornflowerblue', cex=cex_point)
qqconf(length(h2$h2_p_D)+1)
points(sort(-log10(seq(0,1,length.out=(length(h2$h2_p_D)+1))[-1])),
	sort(-log10(h2$h2_p_D)),
	pch=16, col='cornflowerblue', cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)

dev.off()

# Let's take a look at the QQ plots for the various different types of traits

h2_cts <- h2 %>% filter(variable_type == "continuous_irnt")
h2_cat <- h2 %>% filter(variable_type != "continuous_irnt")

pdf("~/Repositories/ldscgxe_dominance_paper/maf_0.01_heritability/heritability_results_summaries/additive_dominance_maf_0.01_qq_plots_cts.pdf", width=9, height=4.5)
par(mfrow=c(1,2))

plot(sort(-log10(seq(0,1,length.out=(length(h2_cts$h2_p_A)+1))[-1])),
	sort(-log10(h2_cts$h2_p_A)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Additive', col='grey20', cex=cex_point)
qqconf(length(h2_cts$h2_p_A)+1)
points(sort(-log10(seq(0,1,length.out=(length(h2_cts$h2_p_A)+1))[-1])),
	sort(-log10(h2_cts$h2_p_A)), pch=16, col='grey20', cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)

plot(sort(-log10(seq(0,1,length.out=(length(h2_cts$h2_p_D)+1))[-1])),
	sort(-log10(h2_cts$h2_p_D)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")),
	ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Dominance', col='grey20', cex=cex_point)
qqconf(length(h2_cts$h2_p_D)+1)
points(sort(-log10(seq(0,1,length.out=(length(h2_cts$h2_p_D)+1))[-1])),
	sort(-log10(h2_cts$h2_p_D)), pch=16, col='grey20', cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)

dev.off()

pdf("~/Repositories/ldscgxe_dominance_paper/maf_0.01_heritability/heritability_results_summaries/additive_dominance_maf_0.01_qq_plots_cat.pdf", width=9, height=4.5)
par(mfrow=c(1,2))
colorBlindGrey8   <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
colours <- colorBlindGrey8[1:4]
plot(sort(-log10(seq(0,1,length.out=(length(h2_cat$h2_p_A)+1))[-1])),
	sort(-log10(h2_cat$h2_p_A)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")), ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Additive', col=colours[1], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 10000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 10000))$h2_p_A)), pch=16, col=colours[2], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 20000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 20000))$h2_p_A)), pch=16, col=colours[3], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 40000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 40000))$h2_p_A)), pch=16, col=colours[4], cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)
legend('topleft',
	legend=c(">3,000 cases", ">10,000 cases", ">20,000 cases", ">40,000 cases"),
	col=colours,
	lty=c(1,1,1,1), bty='n')

plot(sort(-log10(seq(0,1,length.out=(length(h2_cat$h2_p_D)+1))[-1])),
	sort(-log10(h2_cat$h2_p_D)),
	pch=16, xlab=expression(paste("Expected (",-log[10], " p-value)")), ylab=expression(paste("Observed (",-log[10], " p-value)")),
	main='Dominance', col=colours[1], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 10000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 10000))$h2_p_D)), pch=16, col=colours[2], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 20000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 20000))$h2_p_D)), pch=16, col=colours[3], cex=cex_point)
points(sort(-log10(seq(0,1,length.out=(nrow(h2_cat %>% filter(n_cases > 40000))+1))[-1])),
	sort(-log10((h2_cat %>% filter(n_cases > 40000))$h2_p_D)), pch=16, col=colours[4], cex=cex_point)
lines(c(-10,10), c(-10,10), col='red', lwd=line_width)
legend('topleft',
	legend=c(">3,000 cases", ">10,000 cases", ">20,000 cases", ">40,000 cases"),
	col=colours,
	lty=c(1,1,1,1), bty='n')

dev.off()

# I also want to grab a bunch of numbers here for the abstract (disbtribution of h2 on the liability scale etc).