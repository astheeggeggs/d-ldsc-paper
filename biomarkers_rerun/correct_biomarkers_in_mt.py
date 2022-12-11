import sys
import hail as hl

# Finally, combine this with the biomarkers matrix table.
# First, need to combine the dominance biomarker pipelines into a matrix table - see load_results_matrix_table_biomarkers_dominance.py in this folder
mt_dom = hl.read_matrix_table('gs://ukbb-dominance/biomarkers/results/combined/ukb31063.biomarker_gwas_results.both_sexes_rerun.mt')
mt_add = hl.read_matrix_table('gs://ukb31063-mega-gwas/biomarkers/results-matrix-tables/ukb31063.biomarker_gwas_results.both_sexes.mt')

mt_add = mt_add.annotate_entries(AF=mt_add['sum_x'] / (2 * mt_add['n']))
mt_add = mt_add.annotate_entries(
        minor_AF=hl.cond(mt_add['AF'] <= 0.5, mt_add['AF'], 1.0 - mt_add['AF']),
        minor_allele=hl.cond(mt_add['AF'] <= 0.5, mt_add['alleles'][1], mt_add['alleles'][0])
        )
mt_add = mt_add.annotate_entries(low_confidence_variant=mt_add['minor_AF'] < 0.001)

mt_add = mt_add.annotate_entries(additive=hl.struct(
            expected_case_minor_AC = hl.null(hl.tfloat64),
            expected_min_category_minor_AC = hl.null(hl.tfloat64), 
            low_confidence_variant=mt_add['low_confidence_variant'],
            n_complete_samples=mt_add['n'],
            AC=mt_add['sum_x'],
            ytx=mt_add['y_transpose_x'],
            beta=mt_add['beta'],
            se=mt_add['standard_error'],
            tstat=mt_add['t_stat'],
            pval=mt_add['p_value']))
mt = mt_add.select_entries(mt_add.additive)

mt_dom = mt_dom.annotate_entries(dominance=hl.struct(
            AC=mt_dom['AC'],
            ytx=mt_dom['ytx'],
            beta=mt_dom['beta'],
            se=mt_dom['se'],
            tstat=mt_dom['tstat'],
            pval=mt_dom['pval']))
mt_dom = mt_dom.select_entries(mt_dom.dominance)

mt = mt.annotate_entries(dominance = mt_dom[mt.row_key, mt.col_key].dominance)
mt = mt.annotate_cols(PHESANT_transformation = hl.null(hl.tstr),
                      n_cases = hl.null(hl.tint32),
                      n_controls = hl.null(hl.tint32),
                      variable_type = hl.cond(mt.phenotype.contains('irnt'), 'continuous_irnt', 'continuous_raw'),
                      source = 'phesant',
                      notes = mt.description,
                      description = mt.biomarker).rename({'phenotype': 'phenotype_tmp'})
mt = mt.annotate_cols(phenotype = mt.code).key_cols_by('phenotype').drop('phenotype_tmp', 'biomarker', 'code')

mt_existing = hl.read_matrix_table("gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_corrected_jan2020.mt")

temp = mt_existing.key_cols_by()
mt_existing = temp.select_cols(
    **{k: temp[k] for k in mt.col}
    ).key_cols_by(*mt_existing.col_key)

mt = mt_existing.union_cols(mt)

mt.write("gs://ukbb-dominance/additive_dominance_gwas_results.both_sexes_including_nov2022_biomarkers_corrected_jan2020.mt", overwrite = True)

