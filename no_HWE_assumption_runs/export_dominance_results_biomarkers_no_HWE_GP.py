import sys
import hail as hl

# hailctl dataproc submit dp export_dominance_results_biomarkers_no_HWE.py pipeline

pipeline = sys.argv[1]

sex = 'both_sexes'
contig = 'autosomes'

if sex not in set(['both_sexes', 'female', 'male']):
    raise ValueError(f'Invalid sex argument "{sex}" - must be one of {{"both_sexes", "female", "male"}}.')
if contig not in set(['autosomes', 'chrX', 'chrXY']):
    raise ValueError(f'Invalid contig argument "{contig}" - must be one of {{"autosomes", "chrX", "chrXY"}}.')

codes = {
    'albumin': '30600',
    'alkaline_phosphatase': '30610',
    'alanine_aminotransferase': '30620',
    'apoliprotein_A': '30630',
    'apoliprotein_B': '30640',
    'aspartate_aminotransferase': '30650',
    'direct_bilirubin': '30660',
    'urea': '30670',
    'calcium': '30680',
    'cholesterol': '30690',
    'creatinine': '30700',
    'C_reactive_protein': '30710',
    'cystatin_C': '30720',
    'gamma_glutamyltransferase': '30730',
    'glucose': '30740',
    'glycated_haemoglobin': '30750',
    'hdl_cholesterol': '30760',
    'igf_1': '30770',
    'ldl': '30780',
    'lipoprotein_A': '30790',
    'oestradiol': '30800',
    'phosphate': '30810',
    'rheumatoid_factor': '30820',
    'shbg': '30830',
    'total_bilirubin': '30840',
    'testosterone': '30850',
    'total_protein': '30860',
    'triglycerides': '30870',
    'urate': '30880',
    'vitamin_D': '30890',
    'estimated_sample_dilution_factor': '30897'}
    
ht_results = hl.read_table(f'gs://ukbb-dominance/biomarkers/results/no_HWE_assumption/ukb31063.biomarker_gwas_results.{sex}.{contig}.pipeline_{pipeline}_no_HWE_assumption_GP.ht')

ht_results = ht_results.annotate(variant=hl.delimit(hl.array([
    ht_results['locus'].contig,
    hl.str(ht_results['locus'].position),
    ht_results['alleles'][0],
    ht_results['alleles'][1]]), delimiter=':'))
ht_results = ht_results.key_by('variant')
ht_results = ht_results.repartition(116)
ht_results = ht_results.cache()

phenotypes = ht_results['phenotypes'].collect()[0]
for i, phenotype in enumerate(phenotypes):
    variable_type = phenotype.split('_')[-1]
    code = codes[phenotype.replace('_raw', '').replace('_irnt', '')]
    ht_export = ht_results.annotate(
        dominance_AC=ht_results['sum_x'][i],
        dominance_ytx=ht_results['y_transpose_x'][i][0],
        dominance_beta=ht_results['beta'][i][0],
        dominance_se=ht_results['standard_error'][i][0],
        dominance_tstat=ht_results['t_stat'][i][0],
        dominance_pval=ht_results['p_value'][i][0])
    ht_export = ht_export.select(
        'dominance_AC',
        'dominance_ytx',
        'dominance_beta',
        'dominance_se',
        'dominance_tstat',
        'dominance_pval')

    ht_export.export(f'gs://ukb-mega-gwas-results/round2/dominance-biomarkers-tsvs-no-HWE-assumption/{code}_{variable_type}.gwas.imputed_v3.{sex}._no_HWE_assumption_GP.tsv.bgz')
