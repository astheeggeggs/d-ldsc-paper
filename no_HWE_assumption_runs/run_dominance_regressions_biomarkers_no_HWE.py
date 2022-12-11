import sys
import hail as hl

# hailctl dataproc submit dp run_dominance_regressions_biomarkers_no_HWE.py pipeline

hl.init()

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

    if contig == 'autosomes':
        contig_expr = 'chr{1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22}'
    else:
        contig_expr = contig

    mt = hl.import_bgen(
        path=f'gs://fc-7d5088b4-7673-45b5-95c2-17ae00a04183/imputed/ukb_imp_{contig_expr}_v3.bgen',
        sample_file=f'gs://ukb31063/ukb31063.{contig}.sample',
        entry_fields=['GT'],
        variants=ht_variants)

    mt = mt.annotate_cols(
        phenotypes=ht_phenotypes[mt.s],
        covariates=ht_covariates[mt.s]
        )
    
    mt = mt.filter_cols(hl.is_defined(mt.covariates.PC1))

    # Run orthogonalised encoding
    mt = mt.annotate_rows(
        r_GT = hl.agg.mean(mt.GT.is_hom_ref()),
        h_GT = hl.agg.mean(mt.GT.is_het_ref()),
        a_GT = hl.agg.mean(mt.GT.is_hom_var())
    )

    mt = mt.annotate_rows(
        c_D_GT = 1 / hl.sqrt(mt.h_GT * mt.a_GT * mt.r_GT * ( mt.a_GT + mt.r_GT - (mt.a_GT - mt.r_GT) ** 2)),
        ha_GT = mt.h_GT * mt.a_GT,
        ar_GT = mt.a_GT * mt.r_GT,
        hr_GT = mt.h_GT * mt.r_GT,
    )

    mt = mt.annotate_entries(
        x_D_GT = mt.c_D_GT * (
            - mt.ha_GT * mt.GT.is_hom_ref() + 
            2 * mt.ar_GT * mt.GT.is_het_ref() + 
            - mt.hr_GT * mt.GT.is_hom_var()
        )
    )
    
    ht = hl.linear_regression_rows(
        y=[[mt['phenotypes'][y]] for y in phenotypes_to_analyse],
        x=mt.x_D_GT,
        covariates=[1.0, *[mt['covariates'][x] for x in list(mt['covariates'].keys())]],
        pass_through=['varid', 'rsid'])

    ht = ht.annotate_globals(phenotypes = phenotypes_to_analyse)

    if dilution:
        ht.write(f'gs://ukbb-dominance/biomarkers/results/no_HWE_assumption/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}.dilution_factor_no_HWE_assumption.ht',
            overwrite=True)
    else:
        ht.write(f'gs://ukbb-dominance/biomarkers/results/no_HWE_assumption/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}_no_HWE_assumption.ht',
            overwrite=True)

