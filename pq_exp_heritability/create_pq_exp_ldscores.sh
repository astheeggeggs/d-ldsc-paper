#!/bin/bash

# Local and cluster versions of the plink files used to generate the LD-scores.
# scp ~/Repositories/ldscgxe/1000G_EUR_Phase3_plink/1000G.EUR.QC.* dpalmer@login01:/psych/genetics_data/dpalmer/ldscores-pq-exp/1000G_EUR_Phase3_plink/

conda activate d-ldsc
hm3_list="../ldscgxe/1000G_EUR_Phase3_plink/w_hm3.snplist"
cd ~/Repositories/d-ldsc

for chr in {1..22}; do
	bfile="../ldscgxe/1000G_EUR_Phase3_plink/1000G.EUR.QC.$chr"
	for alpha in {-1.0,-0.5,-0.25,0.0,0.25,0.5,1.0}; do
	out="../ldscgxe_dominance_paper/pq_exp_heritability/ldscores/1000G.EUR.QC.alpha.${alpha}.chr${chr}"
	./get_ldscores.py --additive --dominance \
	--ld-wind-cm 1 --extract $hm3_list \
	--bfile $bfile --out $out --pq-exp $alpha
	done
done

# Copy the ldscore files and M_5_50, M, up to the cluster.
# scp ~/Repositories/ldscgxe_dominance_paper/pq_exp_heritability/ldscores/1000G.EUR.QC.alpha* dpalmer@login01:/psych/genetics_data/dpalmer/ldscores-pq-exp/ldscores/
