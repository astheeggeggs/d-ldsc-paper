library(dplyr)
library(ggplot2)
library(data.table)
library(tidyr)
library(lattice)
library(ggrastr)

read_in_sumstats_and_restrict <- function(path_to_results_add, path_to_results_dom, MAF_filter=0.05)
{
  df_add <- fread(paste0('zcat ', path_to_results_add), index='variant', key='variant',
    select=c('variant', 'low_confidence_variant', 'pval', 'minor_AF'))

  df_dom <- fread(paste0('zcat ', path_to_results_dom), index='variant', key='variant',
    select=c('variant', 'dominance_pval'))

  df <- merge(df_add, df_dom) %>% filter(minor_AF > MAF_filter) %>% filter(!low_confidence_variant)

  return(df)
}

make_manhattan_plot <- function(df, df_contigs, plot_title="", include_X=TRUE) {

    if(!include_X) {
      shifted_position <- df_contigs$shifted_position[-23]
      contig_labels <- df_contigs$contig[-23]
    } else {
      shifted_position <- df_contigs$shifted_position
      contig_labels <- gsub(23, 'X', df_contigs$contig)
    }

    p <- ggplot(df) +
    geom_point(size=0.5, color='#4E79A7', aes(x=xvar, y=yvar)) +
    geom_hline(yintercept=-log10(5e-8), color='#F28E2B', linetype='dashed', size=1) +
    scale_x_continuous(breaks=shifted_position, labels=contig_labels) +
    scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
    labs(x='Chromosome', y='-log10(pval)', title=plot_title) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))

    return(p)
}


loglog <- function(x) {
  x[which(x >= 30)] <- 30 + log10(x[which(x>=30)]-30+1)
  return(x)
}

logloginv <- function(x) {
  x[which(x >= 30)] <- 10^(x[which(x>=30)]-30)+30-1
  return(x)
}

loglog2 <- function(x) {
  x[which(x >= 300)] <- 300 + 10*log10(x[which(x>=300)]-300+1)
  return(x)
}

logloginv2 <- function(x) {
  x[which(x >= 300)] <- 10^((x[which(x>=300)] - 300)/10)+300-1
  return(x)
}

loglogtrans <- scales::trans_new('loglog', loglog, logloginv)
loglogtrans2 <- scales::trans_new('loglog2', loglog2, logloginv2)

make_manhattan_plot_coloured <- function(dt, dt_contigs, plot_title="", colour_col, include_X=TRUE, transform=NULL, breaks=NULL) {

    if(!include_X) {
      shifted_position <- dt_contigs$shifted_position[-23]
      contig_labels <- dt_contigs$contig[-23]
    } else {
      shifted_position <- dt_contigs$shifted_position
      contig_labels <- gsub(23, 'X', dt_contigs$contig)
    }

    if(is.null(transform)) {
      p <- ggplot(dt) +
        geom_point_rast(size=0.5, aes(x=xvar, y=yvar, color=factor(colour_col))) +
        geom_hline(yintercept=-log10(5e-8), color='#F28E2B', linetype='dashed', size=1) +
        scale_x_continuous(breaks=shifted_position, labels=contig_labels) +
        scale_y_continuous(breaks=scales::pretty_breaks(n=10)) +
        labs(x='Chromosome', y='-log10(pval)', title=plot_title) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
    } else {
      p <- ggplot(dt) +
        geom_point_rast(size=0.5, aes(x=xvar, y=yvar, color=factor(colour_col))) +
        geom_hline(yintercept=-log10(5e-8), color='darkgray', size=0.5) +
        scale_x_continuous(breaks=shifted_position, labels=contig_labels) +
        scale_y_continuous(trans=transform, breaks=breaks) +
        labs(x='Chromosome', y='-log10(pval)', title=plot_title) +
        theme_bw() +
        theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank())
    }

    return(p)
}

make_ukbb_summary_manhattan_plot <- function(path_to_results, tstat_col, colour_col, plot_title="", plot_name="out",
  MAF_filter=0.05, p_HWE_filter=1e-6, log_p_val_threshold=3, transform=NULL, breaks=NULL)
{
  buffer <- 1e8
  contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
  dt_contigs <- data.frame(contig=contigs_,
                          start=c(1, 249250623, 492449997, 690472428,  881626705, 1062541966, 1233657034, 1392795698,
                                  1539159721,  1680373153, 1815907901, 1950914418, 2084766314, 2199936193, 2307285734,
                                  2409817127,  2500171881, 2581367092, 2659444341, 2718573325, 2781598846, 2829728742,
                                  2881033309),
                          end=c( 249250622,  492449996,  690472427,  881626704, 1062541965, 1233657033, 1392795697,
                                1539159720, 1680373152, 1815907900, 1950914417, 2084766313, 2199936192, 2307285733,
                                2409817126, 2500171880, 2581367091, 2659444340, 2718573324, 2781598845, 2829728741,
                                2881033308, 3036303869)) %>%
  mutate(middle=floor(start + (end - start)/2),
         length=(end - start)) %>%
  mutate(shifted_position=(middle + (contig - 1)*buffer))

  if (length(grep("gs://", path_to_results)) > 0) {
    read_in <- paste0("gsutil cat ", path_to_results, " | gzcat")
  } else {
    read_in <- paste0("gzcat ", path_to_results)
  }

  dt <- fread(read_in, index='variant', key='variant',
    select=c('variant', 'p_hwe', 'minor_AF', tstat_col, colour_col)) %>% filter(minor_AF > MAF_filter, p_hwe > p_HWE_filter)

  dt <- dt %>% mutate(negative_log_p = -((log(2) + pnorm(abs(dt[[tstat_col]]), log.p=TRUE, lower.tail=FALSE)) / log(10)),
    colour_col = dt[[colour_col]])
  
  dt <- dt %>% filter(negative_log_p > log_p_val_threshold) %>%
  separate(variant, c('contig', 'pos', 'ref', 'alt'), sep=':', remove=FALSE) %>% filter(contig!='X') %>%
  mutate(chrom=gsub('X', '23', contig), pos=as.integer(pos)) %>%
  mutate(xvar=dt_contigs[gsub('X', '23', contig), 'start'] + pos + (as.integer(gsub('X', '23', contig))-1)*buffer,
         yvar=negative_log_p) %>%
  select(xvar, yvar, colour_col)

  p <- make_manhattan_plot_coloured(dt, dt_contigs, plot_title=plot_title, include_X=FALSE, transform=transform, breaks=breaks)
  
  ggsave(paste0(plot_name,'.pdf'), plot=p, width=50, height=9, units='cm')
}

make_ukbb_manhattan_plot <- function(path_to_results_add, path_to_results_dom, plot_title="", plot_name="out", MAF_filter=0.05)
{
  buffer <- 1e8
  contigs_ <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)
  df_contigs <- data.frame(contig=contigs_,
                          start=c(1, 249250623, 492449997, 690472428, 881626705, 1062541966, 1233657034, 1392795698,
                                  1539159721, 1680373153, 1815907901, 1950914418, 2084766314, 2199936193, 2307285734, 2409817127,
                                  2500171881, 2581367092, 2659444341,2718573325, 2781598846, 2829728742, 2881033309),
                          end=c(249250622, 492449996, 690472427, 881626704, 1062541965, 1233657033, 1392795697, 1539159720,
                                1680373152, 1815907900, 1950914417, 2084766313, 2199936192, 2307285733, 2409817126, 2500171880,
                                2581367091, 2659444340, 2718573324, 2781598845, 2829728741, 2881033308, 3036303869)) %>%
  mutate(middle=floor(start + (end - start)/2),
         length=(end - start)) %>%
  mutate(shifted_position=(middle + (contig - 1)*buffer))

  df <- read_in_sumstats_and_restrict(path_to_results_add, path_to_results_dom, MAF_filter)
  
  df_add <- df %>% filter(pval <= 0.01) %>%
  separate(variant, c('contig', 'pos', 'ref', 'alt'), sep=':', remove=FALSE) %>%
  mutate(chrom=gsub('X', '23', contig), pos=as.integer(pos)) %>%
  mutate(xvar=df_contigs[gsub('X', '23', contig), 'start'] + pos + (as.integer(gsub('X', '23', contig))-1)*buffer,
         yvar=-log10(pval)) %>%
  select(xvar, yvar)

  df_dom <- df %>% filter(dominance_pval <= 0.01) %>%
  separate(variant, c('contig', 'pos', 'ref', 'alt'), sep=':', remove=FALSE) %>%
  filter(contig != 'X') %>%
  mutate(chrom=gsub('X', '23', contig), pos=as.integer(pos)) %>%
  mutate(xvar=df_contigs[gsub('X', '23', contig), 'start'] + pos + (as.integer(gsub('X', '23', contig))-1)*buffer,
         yvar=-log10(dominance_pval)) %>%
  select(xvar, yvar)

  p_add <- make_manhattan_plot(df_add, df_contigs, plot_title=plot_title)
  p_dom <- make_manhattan_plot(df_dom, df_contigs, plot_title=plot_title, include_X=FALSE)
  
  ggsave(paste0(plot_name,'_add.pdf'), plot=p_add, width=25, height=9, units='cm')
  ggsave(paste0(plot_name,'_dom.pdf'), plot=p_dom, width=25, height=9, units='cm')
}

make_ukbb_manhattan_plots <- function(path_to_pairs, output_path, MAF_filter=0.05) {
  # Should also include the name of the phenotype and N_cases/N_controls for each of the plots, or N.
  current_pairs <- fread(path_to_pairs)
  current_phenos <- gsub("^.*\\/", "", current_pairs$X1)
  current_phenos <- sapply(strsplit(current_phenos, split='\\.'), `[[`, 1)
  for(i in 1:nrow(current_pairs)) {
    make_ukbb_manhattan_plot(current_pairs[i,1], current_pairs[i,2],
      plot_title=current_phenos[i], plot_name=paste0(output_path,'/', current_phenos[i]),
      MAF_filter=MAF_filter)
  }
}

make_ukbb_qq_plots <- function(path_to_pairs, output_path, filter=TRUE) {
  # Should also include the name of the phenotype and N_cases/N_controls for each of the plots, or N.
  current_pairs <- fread(path_to_pairs)
  current_phenos <- gsub("^.*\\/", "", current_pairs$X1)
  current_phenos <- sapply(strsplit(current_phenos, split='\\.'), `[[`, 1)
  for(i in 1:nrow(current_pairs)) {
    create_qq(current_pairs[i,1], current_pairs[i,2],
      plot_title=current_phenos[i], file_out=paste0(output_path,'/', current_phenos[i]),
      filter=filter)
  }
}

get_files <- function(path_to_files) {
  sumstats <- system(paste0('ls ', path_to_files), intern=TRUE)
  sumstats <- sumstats[grep(".tsv.bgz$", sumstats)]
  return(sumstats)
}
# Want to get a collection of 1000 pairs, and write this to disk.
# Additive results:
# /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/{both_sexes,males,females}/
# Dominance results:
# /psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/{both_sexes,males,females}/

get_list_of_pairs <- function(path_to_results_add, path_to_results_dom,
  results_add_tail=".gwas.imputed_v3.male.tsv.bgz", results_dom_tail=".dominance.gwas.imputed_v3.male.tsv.bgz",
  max_num_pairs=333, out='males') {
  sumstats_add <- get_files(path_to_results_add)
  sumstats_dom <- get_files(path_to_results_dom)
  pheno_add <- sapply(strsplit(sumstats_add, split='\\.'), `[[`, 1)
  pheno_dom <- sapply(strsplit(sumstats_dom, split='\\.'), `[[`, 1)
  phenos <- intersect(pheno_add, pheno_dom)
  file_pairs <- cbind(paste0(path_to_results_add, phenos, results_add_tail),
                      paste0(path_to_results_dom, phenos, results_dom_tail))
  # Export collection of 333 pairs
  num_per_run <- ceiling(length(phenos) / max_num_pairs)
  num_pairs <- ceiling(length(phenos) / num_per_run)
  for(i in 1:(num_pairs-1)) {
    fwrite(data.frame(file_pairs[((i-1)*num_per_run+1):(i*num_per_run),]), file=paste0(out, ".", i, ".phenos.txt"))
  }
  fwrite(data.frame(file_pairs[((num_pairs-1)*num_per_run+1):nrow(file_pairs),]), file=paste0(out, ".", num_pairs, ".phenos.txt"))
}

create_qq <- function(sumstats_file_add, sumstats_file_dom, plot_title, file_out, filter=TRUE, write_pdf=TRUE)
{
  sumstats_add <- fread(paste0('zcat ', sumstats_file_add),
    index='variant', key='variant', select = c("variant", "pval", "low_confidence_variant", "minor_AF"))
  sumstats_dom <- fread(paste0('zcat ', sumstats_file_dom),
    index='variant', key='variant', select = c("variant", "dominance_pval"))
  sumstats <- merge(sumstats_add, sumstats_dom, by="variant")
  
  if(filter == TRUE)
    sumstats <- sumstats %>% filter(low_confidence_variant == FALSE) %>% filter(minor_AF > 0.05)
  
  # Sample a subset.

  if(write_pdf==TRUE) {
    pdf(file=file_out)
  }
  create_qq_plot(sumstats$pval, paste(plot_title, "additive"))
  create_qq_plot(sumstats$dominance_pval[-grep('X', sumstats$variant)], paste(plot_title, "dominance"))
  if(write_pdf==TRUE) {
    dev.off()
  }
}

create_qq_plot <- function(pvalues, title, should.thin=TRUE, thin.obs.places=2, thin.exp.places=2, 
  xlab=expression(paste("Expected (", -log[10], " p-value)")),
  ylab=expression(paste("Observed (", -log[10], " p-value)")), 
  draw.conf=TRUE, conf.points=1000, conf.col="lightgray", conf.alpha=.05,
  already.transformed=FALSE, pch=20, aspect="iso",
  par.settings=list(superpose.symbol=list(pch=pch)), ...)
{
  
  if(any(is.na(pvalues)))
    pvalues <- pvalues[-is.na(pvalues)]

  if(any(pvalues == 0))
    pvalues[pvalues==0] <- .Machine$double.xmin

  n <- length(pvalues)+1
  exp.x <- -log10((rank(pvalues, ties.method="first") - 0.5)/n)
  pvalues <- -log10(pvalues)

  # This is a helper function to draw the confidence interval
  qqconf <- function(n, conf.points=1000, conf.col="gray", conf.alpha=.05)
  {
    conf.points = min(conf.points, n-1)
    mpts <- matrix(nrow=conf.points*2, ncol=2)
    
    seq_points <- seq(from=1, to=conf.points)
    mpts[seq_points,1] <- -log10((seq_points-0.5)/n)
    mpts[seq_points,2] <- -log10(qbeta(1-conf.alpha/2, seq_points, n-seq_points))
    mpts[conf.points*2+1-seq_points,1] <- -log10((seq_points-0.5)/n)
    mpts[conf.points*2+1-seq_points,2] <- -log10(qbeta(conf.alpha/2, seq_points, n-seq_points))

    polygon(x=mpts[,1], y=mpts[,2], border=FALSE, col='grey')
  }

  # Reduce number of points to plot
  thin <- unique(data.frame(pvalues = round(pvalues, thin.obs.places), exp.x = round(exp.x, thin.exp.places)))
  pvalues <- thin$pvalues
  exp.x <- thin$exp.x

  # Draw the plot
  plot(c(0,max(exp.x, na.rm=TRUE)), c(0, max(pvalues, na.rm=TRUE)), pch=pch, xlab=xlab, ylab=ylab, main=title)
  qqconf(n)
  points(exp.x, pvalues, pch=pch)
  lines(c(0, max(exp.x, na.rm=TRUE)), c(0, max(exp.x, na.rm=TRUE)), col='#F28E2B')

}
