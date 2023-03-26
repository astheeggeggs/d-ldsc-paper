#\!/bin/bash
#$ -N md5sum
#$ -cwd
#$ -w e
#$ -j y
#$ -b n
#$ -l h_vmem=16g
#$ -l h_rt=72:00:00
#$ -l os=RedHat7

source /broad/software/scripts/useuse
use UGER

md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/*tsv.bgz > both_sexes_md5sum_cluster.tsv
md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/biomarkers/*tsv.bgz >> both_sexes_md5sum_cluster.tsv

md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/males/*tsv.bgz > male_md5sum_cluster.tsv
md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/males/biomarkers/*tsv.bgz >> male_md5sum_cluster.tsv

md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/females/*tsv.bgz > female_md5sum_cluster.tsv
md5sum /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/females/biomarkers/*tsv.bgz >> female_md5sum_cluster.tsv
