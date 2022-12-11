import hail as hl

mt = hl.read_matrix_table('gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_including_nov2022_biomarkers_corrected_jan2020.mt')

# Read in and annotate with the data we have generated in `ldscgxe_dominance_paper/generate_raymonds_primary_both_sex_phenotypes.r`.
# gsutil cp additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv gs://ukbb-dominance/

ht_pheno = hl.import_table('gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.phenotypes_and_boolean_filters.tsv', impute=True)
ht_pheno = ht_pheno.key_by('phenotype')
# Check the count.

mt = mt.annotate_cols(
  phenotype_qc = ht_pheno[mt.phenotype].raymond_both_sexes_qc,
  dominance_phenotype_qc = ht_pheno[mt.phenotype].dominance_qc,
  dominance_phenotype_qc_no_ordinal = ht_pheno[mt.phenotype].dominance_qc_no_ordinal)

mt = mt.filter_cols(mt.dominance_phenotype_qc == True)
mt_no_ord = mt.filter_cols(mt.dominance_phenotype_qc_no_ordinal == True)

# Sanity checks - these should do nothing.
# Remove phenotypes that don't have sufficient data
mt = mt.filter_cols((hl.is_defined(mt.n_cases) & (mt.n_cases > 3000) & (mt.n_controls > 3000)) | 
                    (~hl.is_defined(mt.n_cases) & (mt.n_non_missing > 50000)))
# Remove raw cts phenotypes
mt = mt.filter_cols(~mt.phenotype.contains('raw'))
mt_no_ord = mt.filter_cols(~mt.variable_type.contains('ordinal'))

mt = mt.annotate_rows(extreme_dominance_tstat = hl.agg.filter(
    ~hl.is_nan(mt.dominance.tstat), hl.agg.max(hl.abs(mt.dominance.tstat))),
                     extreme_additive_tstat = hl.agg.filter(
    ~hl.is_nan(mt.additive.tstat), hl.agg.max(hl.abs(mt.additive.tstat))))
mt.count()

mt = mt.annotate_rows(extreme_phenotype_dominance = hl.agg.filter(hl.abs(mt.dominance.tstat) == mt.extreme_dominance_tstat, hl.agg.take(mt.phenotype, 1)))
mt = mt.annotate_rows(extreme_phenotype_additive = hl.agg.filter(hl.abs(mt.additive.tstat) == mt.extreme_additive_tstat, hl.agg.take(mt.phenotype, 1)))
# A few don't match exactly due to floating point precision.
mt = mt.filter_rows((mt.extreme_phenotype_dominance.length() > 0) & (mt.extreme_phenotype_additive.length() > 0))
mt.count()

mt = mt.annotate_rows(extreme_phenotype_additive = mt.extreme_phenotype_additive[0],
                      extreme_phenotype_dominance = mt.extreme_phenotype_dominance[0])

mt = mt.annotate_rows(extreme_additive_struct = hl.agg.filter(hl.abs(mt.additive.tstat) == mt.extreme_additive_tstat, hl.agg.take(mt.additive, 1)[0]))
mt = mt.annotate_rows(extreme_dominance_struct = hl.agg.filter(hl.abs(mt.dominance.tstat) == mt.extreme_dominance_tstat, hl.agg.take(mt.dominance, 1)[0]))

field_names_A = [x for x in mt.extreme_additive_struct]
field_names_D = [x for x in mt.extreme_dominance_struct]

mt = mt.annotate_rows(**{"extreme_additive_" + field_name : mt.extreme_additive_struct[field_name] for field_name in field_names_A},
                      **{"extreme_dominance_" + field_name : mt.extreme_dominance_struct[field_name] for field_name in field_names_D})

mt = mt.drop(mt.extreme_additive_struct, mt.extreme_dominance_struct)
mt.rows().export('gs://ukbb-dominance/aggregate_manhattan_rerun_biomarkers.tsv.bgz')
mt.filter_rows(mt.minor_AF > 0.05).rows().export('gs://ukbb-dominance/aggregate_manhattan_maf_0_05_rerun_biomarkers.tsv.bgz')

# Do the same excluding ordinal variables.
mt_no_ord = mt_no_ord.annotate_rows(extreme_dominance_tstat = hl.agg.filter(
    ~hl.is_nan(mt_no_ord.dominance.tstat), hl.agg.max(hl.abs(mt_no_ord.dominance.tstat))),
                     extreme_additive_tstat = hl.agg.filter(
    ~hl.is_nan(mt_no_ord.additive.tstat), hl.agg.max(hl.abs(mt_no_ord.additive.tstat))))
mt_no_ord.count()

mt_no_ord = mt_no_ord.annotate_rows(extreme_phenotype_dominance = hl.agg.filter(hl.abs(mt_no_ord.dominance.tstat) == mt_no_ord.extreme_dominance_tstat, hl.agg.take(mt_no_ord.phenotype, 1)))
mt_no_ord = mt_no_ord.annotate_rows(extreme_phenotype_additive = hl.agg.filter(hl.abs(mt_no_ord.additive.tstat) == mt_no_ord.extreme_additive_tstat, hl.agg.take(mt_no_ord.phenotype, 1)))
# A few don't match exactly due to floating point precision.
mt_no_ord = mt_no_ord.filter_rows((mt_no_ord.extreme_phenotype_dominance.length() > 0) & (mt_no_ord.extreme_phenotype_additive.length() > 0))
mt_no_ord.count()

mt_no_ord = mt_no_ord.annotate_rows(extreme_phenotype_additive = mt_no_ord.extreme_phenotype_additive[0],
                                    extreme_phenotype_dominance = mt_no_ord.extreme_phenotype_dominance[0])

mt_no_ord = mt_no_ord.annotate_rows(extreme_additive_struct = hl.agg.filter(hl.abs(mt_no_ord.additive.tstat) == mt_no_ord.extreme_additive_tstat, hl.agg.take(mt_no_ord.additive, 1)[0]))
mt_no_ord = mt_no_ord.annotate_rows(extreme_dominance_struct = hl.agg.filter(hl.abs(mt_no_ord.dominance.tstat) == mt_no_ord.extreme_dominance_tstat, hl.agg.take(mt_no_ord.dominance, 1)[0]))

field_names_A = [x for x in mt_no_ord.extreme_additive_struct]
field_names_D = [x for x in mt_no_ord.extreme_dominance_struct]

mt_no_ord = mt_no_ord.annotate_rows(**{"extreme_additive_" + field_name : mt_no_ord.extreme_additive_struct[field_name] for field_name in field_names_A},
                                    **{"extreme_dominance_" + field_name : mt_no_ord.extreme_dominance_struct[field_name] for field_name in field_names_D})

mt_no_ord = mt_no_ord.drop(mt_no_ord.extreme_additive_struct, mt_no_ord.extreme_dominance_struct)
mt_no_ord.rows().export('gs://ukbb-dominance/aggregate_manhattan_no_ord_rerun_biomarkers.tsv.bgz')
mt_no_ord.filter_rows(mt_no_ord.minor_AF > 0.05).rows().export('gs://ukbb-dominance/aggregate_manhattan_no_ord_maf_0_05_rerun_biomarkers.tsv.bgz')

# MAF 0.01 cutoff
ht_no_ord = hl.import_table('gs://ukbb-dominance/aggregate_manhattan_no_ord_rerun_biomarkers.tsv.bgz')
ht_no_ord.filter(ht_no_ord.minor_AF > 0.01).export('gs://ukbb-dominance/aggregate_manhattan_no_ord_maf_0_01_rerun_biomarkers.tsv.bgz')

# Below to include inbreeding coefficients (for artefact checking of dominance signal).

ht = hl.import_table('gs://ukbb-dominance/aggregate_manhattan.tsv.bgz', impute=True)
ht = ht.annotate(**hl.parse_variant(ht.variant))
ht = ht.key_by(ht.locus, ht.alleles)

bgen = hl.import_bgen('gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}_v3.bgen',
                      sample_file='gs://ukb31063/ukb31063.autosomes.sample',
                      entry_fields=['GT'])
samples = hl.read_table("gs://ukb31063/ukb31063.neale_gwas_samples.both_sexes.ht")
variants = hl.read_table("gs://ukb31063/ukb31063.neale_gwas_variants.ht")

# Filter to the rows and columns in the GWASes that we ran.
bgen = bgen.filter_cols(hl.is_defined(samples[bgen.s]))
bgen = bgen.filter_rows(hl.is_defined(variants[bgen.row_key]))

# Filter to just the variants that are in the aggregrate GWAS file.
# Annotate with all the information from the overlay of GWAS'.
bgen = bgen.filter_rows(hl.is_defined(ht[bgen.row_key]))
bgen = bgen.annotate_rows(**ht[bgen.row_key])
bgen = bgen.annotate_rows(IB = hl.agg.inbreeding(bgen.GT, bgen.AF))

bgen_rows = bgen.rows()
bgen_rows.select(
    'rsid', 'varid', 'variant', 'minor_allele', 'minor_AF',
    'consequence', 'consequence_category', 'info', 'call_rate',
    'alt_AC', 'AF', 'p_hwe', 'n_not_called', 
    'n_hom_ref', 'n_het', 'n_hom_var', 'n_non_ref',
    'r_heterozygosity', 'r_het_hom_var', 'r_expected_het_frequency',
    'extreme_dominance_tstat', 'extreme_additive_tstat',
    'extreme_dominance_beta', 'extreme_additive_beta',
    'extreme_phenotype_dominance', 'extreme_phenotype_additive',
    **bgen_rows.IB).export("gs://ukbb-dominance/aggregate_manhattan_inbreeding.tsv.bgz")


