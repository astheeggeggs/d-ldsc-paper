import sys
import hail as hl

# hailctl dataproc submit dp run_dominance_regressions_biomarkers.py pipeline

pipeline = sys.argv[1]

sex = 'both_sexes'
contig = 'autosomes'
dilution = False

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')
if contig not in set(['autosomes', 'chrX', 'chrXY']):
    raise ValueError(f'Invalid contig argument "{contig}" - must be one of {{"autosomes", "chrX", "chrXY"}}.')

print(f'starting pipeline {pipeline}')

ht_phenotypes = hl.read_table(f'gs://ukb31063-mega-gwas/biomarkers/pipelines/ukb31063.biomarkers_gwas.{sex}.pipeline_{pipeline}.ht')
    
phenotypes = list(ht_phenotypes[ht_phenotypes.s].keys())

phenotypes_to_analyse = []
for phenotype in phenotypes:
    if phenotype.split('_')[-1] == "irnt":
        phenotypes_to_analyse.append(phenotype)
    
if len(phenotypes_to_analyse) > 0:

    print("analysing")
    print(phenotypes_to_analyse)

    ht_covariates = hl.read_table(f'gs://ukb31063/ukb31063.neale_gwas_covariates.{sex}.ht')
    ht_variants = hl.read_table('gs://ukb31063/ukb31063.neale_gwas_variants.ht')

    if dilution:
        ht = hl.read_table(f'gs://ukb31063/ukb31063.biomarkers_gwas.{sex}.ht')
        ht = ht.select('estimated_sample_dilution_factor_raw')
        ht_covariates = ht_covariates.annotate(
            estimated_sample_dilution_factor=ht[ht_covariates.s]['estimated_sample_dilution_factor_raw'])

    if contig == 'autosomes':
        contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
    else:
        contig_expr = contig

    mt = hl.import_bgen(
        path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
        sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
        entry_fields=['GP'],
        variants=ht_variants)

    mt = mt.annotate_cols(
        phenotypes=ht_phenotypes[mt.s],
        covariates=ht_covariates[mt.s]
        )

    mt = mt.filter_cols((hl.is_defined(mt.covariates.PC1)))

    mt = mt.annotate_rows(p = (hl.agg.sum(mt.GP[1] + 2 * mt.GP[2]) / hl.agg.count_where(hl.is_defined(mt.GP[1] + 2 * mt.GP[2]))) / 2.0)
    mt = mt.annotate_rows(c = mt.p / (mt.p - 1.0))
    mt = mt.annotate_rows(one_over_c = 1.0 / mt.c)
    mt = mt.select_entries(x_D = mt.GP[0] * mt.c + mt.GP[1] + mt.GP[2] * mt.one_over_c)
     
    ht = hl.linear_regression_rows(
        y=[[mt['phenotypes'][y]] for y in phenotypes_to_analyse],
        x=mt.x_D,
        covariates=[1.0, *[mt['covariates'][x] for x in list(mt['covariates'].keys())]],
        pass_through=['varid', 'rsid'])

    ht = ht.annotate_globals(phenotypes=phenotypes_to_analyse)

    if dilution:
        ht.write(f'gs://ukbb-dominance/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}.dilution_factor_rerun.ht',
            overwrite=True)
    else:
        ht.write(f'gs://ukbb-dominance/biomarkers/results/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}_rerun.ht',
            overwrite=True)

