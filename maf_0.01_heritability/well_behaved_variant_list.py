#! /usr/bin/env python

#############
# Based on https://github.com/Nealelab/UK_Biobank_GWAS/blob/master/0.2/create_ldsc_hm3_table.py

# Gets set of SNPs from UKB GWAS to export for ldsc
# - UKB MAF > 1%
# - UKB INFO > 0.9
# - autosomal SNP
# - biallelic (in UKB MAF/INFO filtered set)
# - not in MHC

# this version using Hail 0.2 
# and variant QC from round2 GWAS

# Setup

# re-run everything?
force_all = True

# count variants at each step?
verbose = True

# Input

# UKBB variant QC in the full sample (for INFO scores)
GWAS_mfi = 'gs://ukb31063/ukb31063.neale_gwas_variants.imputed_v3.mfi.ht'

# UKBB variant QC in the GWAS samples (for MAF)
# Note: using saved QC metrics that match values before last round of sample withdrawals
GWAS_qc = 'gs://ukb-mega-gwas-results/round2/annotations/variants.tsv.bgz'

# SNPs passing QC to be included in the UKBB GWAS
GWAS_snps = 'gs://ukb31063/ukb31063.neale_gwas_variants.ht'

### Output

# where to save files
BUCKET = 'gs://ukbb-dominance/ukb_maf_auto_biallelic_maf_0.01_nomhc_snplist/'
stem = 'b37.auto_bi_af_0.01.ukbb_gwas_qcpos.no_mhc'

print("Preparing packages, etc...")

# Filter to biallelic SNPs in UKBB that have MAF > 0.01 and no MHC, and let this intersect with the EUR LD-scores.
# Increase the INFO score to 0.9

# load packages
import hail as hl

# Create filtered keytable with autosomal, biallelic HM3 snps
# with MAF > .01, INFO > 0.9 and passing QC in UKBB GWAS analysis set
    
print("Creating Hail table of SNPs passing UKBB GWAS QC...")

# get list of SNPs to be used in GWAS
# filter here: autosomes, no indels, no MHC (filter early for efficiency)
ukb_snps = hl.read_table(GWAS_snps).key_by('locus','alleles').repartition(500, shuffle=False)

if verbose:
    print("\nCount 1: " + str(ukb_snps.count()) + '\n')

ukb_snps = ukb_snps.filter(
                hl.is_snp(ukb_snps.alleles[0], ukb_snps.alleles[1]) &
                (~(ukb_snps.locus.contig == 'X')) &
                (~((ukb_snps.locus.contig == '6') & (ukb_snps.locus.position > 25000000) & (ukb_snps.locus.position < 34000000)))
            )
if verbose:
    print("\nCount 2: " + str(ukb_snps.count()) + '\n')

# merge in, filter on MAF from the UKBB GWAS sample
ukb_qc = hl.import_table(GWAS_qc)
ukb_qc = ukb_qc.annotate(vstruct = hl.parse_variant(ukb_qc.variant))
ukb_qc = ukb_qc.annotate(locus = ukb_qc.vstruct.locus, alleles = ukb_qc.vstruct.alleles).key_by('locus','alleles')
ukb_qc2 = ukb_snps.join(ukb_qc.select(ukb_qc.minor_AF))
if verbose:
    print("\nCount 3: " + str(ukb_qc2.count()) + '\n')

ukb_qc2 = ukb_qc2.filter(
                (hl.float(ukb_qc2.minor_AF) > 0.01) & 
                (hl.float(ukb_qc2.minor_AF) < 0.99)
            )
if verbose:
    print("\nCount 4: " + str(ukb_qc2.count()) + '\n')

# merge in rsid, info (from full UKB sample)
# and filter to info > 0.9
ukb_mfi = hl.read_table(GWAS_mfi).key_by('locus','alleles').repartition(500, shuffle=False)
ukb_qc3 = ukb_qc2.join(ukb_mfi.select('varid','variant','rsid','info'))
if verbose:
    print("\nCount 5: " + str(ukb_qc3.count()) + '\n')

ukb_qc3 = ukb_qc3.filter(ukb_qc3.info > 0.9)

if verbose:
    print("\nCount 6: " + str(ukb_qc3.count()) + '\n')

# drop multi-allelic sites
loc_count = ukb_qc3.group_by(ukb_qc3.locus).aggregate(nloc=hl.agg.count())
loc_count = loc_count.filter(loc_count.nloc==1)

if verbose:
    print("\nCount 7: " + str(loc_count.count()) + '\n')

ukb_qc3 = ukb_qc3.key_by('locus').join(loc_count).drop('nloc')
if verbose:
    print("\nCount 8: " + str(ukb_qc3.count()) + '\n')

ukb_qc3.write(BUCKET + stem + '.ht', overwrite=True)

# format tsv version
ukb_qc3 = hl.read_table(BUCKET + stem + '.ht')
ukb_qc3 = ukb_qc3.annotate(A1=ukb_qc3.alleles[0],A2=ukb_qc3.alleles[1])
ukb_qc3 = ukb_qc3.key_by('variant')
ukb_qc3 = ukb_qc3.select(ukb_qc3.locus, ukb_qc3.A1, ukb_qc3.A2, ukb_qc3.rsid, ukb_qc3.varid, UKB_GWAS_MAF=ukb_qc3.minor_AF, UKB_all_INFO=ukb_qc3.info)

# save
ukb_qc3.export(BUCKET + stem + '.tsv.bgz')
