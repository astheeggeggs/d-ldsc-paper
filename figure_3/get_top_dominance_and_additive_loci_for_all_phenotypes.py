import hail as hl

main_gs_folder = 'gs://ukbb-dominance'
file_root = 'additive_dominance_gwas_results.both_sexes_including_nov2022_biomarkers_corrected_jan2020'
file_root_tsv = 'additive.dominance.gwas.imputed_v3.both_sexes_including_nov2022_biomarkers_corrected_jan2020_top_hits_' 
# Fig. 3.
# mt = hl.read_matrix_table(f'{main_gs_folder}/additive_dominance_gwas_results.both_sexes_including_biomarkers_corrected_jan2020.mt')
mt = hl.read_matrix_table(f'{main_gs_folder}/{file_root}.mt')

# Name is incorrect, no MAF 0.05 cutoff!
v = hl.import_table(f'{main_gs_folder}/all_mafge0.05.hg19_multianno.csv', delimiter=',', quote = '"', missing='.', impute=True)
v = v.annotate(locus = hl.locus(contig = v.Chr, pos = v.Start), alleles = [v.Ref, v.Alt])
v = v.key_by(v.locus, v.alleles)

v = v.select(v.cytoBand)

mt = mt.annotate_rows(annovar = v[mt.row_key])
mt = mt.annotate_entries(MAF = (mt.additive.AC / (2 * mt.additive.n_complete_samples)))
mt = mt.annotate_entries(MAF = hl.case().when(mt.MAF > 0.5, 1 - mt.MAF)
													.default(mt.MAF))
mt = mt.filter_entries(mt.additive.low_confidence_variant == False)
mt = mt.filter_rows(mt.locus.in_autosome())
mt = mt.checkpoint(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband_checkpoint.mt', overwrite=True)

for MAF_query in [0.01, 0.05]:
    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband_checkpoint.mt')
    if str(MAF_query) == '0.01':
        mt = mt.filter_rows(mt.p_hwe > 1e-10)
    else:
        mt = mt.filter_rows(mt.p_hwe > 1e-6)
    mt.group_rows_by(mt.annovar.cytoBand).aggregate(
        additive = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.additive, 1, ordering = mt.additive.pval)),
        dominance = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.dominance, 1, ordering = mt.dominance.pval))).write(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband.mt', overwrite=True)

    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband.mt')

    mt_final = mt.annotate_cols(
        additive_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].beta, 5, ordering = mt.additive[0].pval)),
        dominance_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].beta, 5, ordering = mt.dominance[0].pval)),
        additive_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].pval, 5, ordering = mt.additive[0].pval)),
        dominance_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].pval, 5, ordering = mt.dominance[0].pval)), 
        AC = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].AC, 5, ordering = mt.additive[0].pval)),
        n_complete_samples = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].n_complete_samples, 5, ordering = mt.additive[0].pval)) 
    )  

    mt_final = mt_final.annotate_cols(MAF = (mt_final.AC / (2 * mt_final.n_complete_samples)))
    mt_final = mt_final.annotate_cols(additive_hits_check = hl.map(
        lambda i: hl.sqrt(2 * mt_final.MAF[i] * (1.0 - mt_final.MAF[i])),
        [0,1,2,3,4]))
    mt_final = mt_final.annotate_cols(additive_hits_rescaled = mt_final.additive_hits / mt_final.additive_hits_check)

    cols = mt_final.cols()
    cols = cols.annotate(**{f'beta_A_{i}': cols.additive_hits_rescaled[i] for i in range(5)})
    cols = cols.annotate(**{f'beta_D_{i}': cols.dominance_hits[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_A_{i}': cols.additive_pval[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_D_{i}': cols.dominance_pval[i] for i in range(5)})

    cols.drop('additive_hits', 'dominance_hits', 'additive_pval', 'dominance_pval').export(f'{main_gs_folder}/{file_root_tsv}no_X_maf_{MAF_query}.tsv.bgz')

# April 2021.
# Next, do the same but filter to additive hits and consider the dominance hits at these same locations.
# Fig. 3 - alternate version, for supplement.
mt = hl.read_matrix_table(f'{main_gs_folder}/{file_root}.mt')

v = hl.import_table(f'{main_gs_folder}/all_mafge0.05.hg19_multianno.csv', delimiter=',', quote = '"', missing='.', impute=True)
v = v.annotate(locus = hl.locus(contig = v.Chr, pos = v.Start), alleles = [v.Ref, v.Alt])
v = v.key_by(v.locus, v.alleles)

v = v.select(v.cytoBand)

mt = mt.annotate_rows(annovar = v[mt.row_key])
mt = mt.annotate_entries(MAF = (mt.additive.AC / (2 * mt.additive.n_complete_samples)))
mt = mt.annotate_entries(MAF = hl.case().when(mt.MAF > 0.5, 1 - mt.MAF)
                                                    .default(mt.MAF))
mt = mt.filter_entries(mt.additive.low_confidence_variant == False)
mt = mt.filter_rows(mt.locus.in_autosome())
mt = mt.checkpoint(f'{main_gs_folder}/group_by_cytoband/{file_root}_additive_hits_by_cytoband_checkpoint.mt', overwrite=True)

for MAF_query in [0.01, 0.05]:
    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband_checkpoint.mt')
    if str(MAF_query) == '0.01':
        mt = mt.filter_rows(mt.p_hwe > 1e-10)
    else:
        mt = mt.filter_rows(mt.p_hwe > 1e-6)
    mt.group_rows_by(mt.annovar.cytoBand).aggregate(
        additive = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.additive, 1, ordering = mt.additive.pval)),
        dominance = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.dominance, 1, ordering = mt.additive.pval))).write(f'{main_gs_folder}/group_by_cytoband/{file_root}_additive_hits_by_cytoband.mt', overwrite=True)

    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_additive_hits_by_cytoband.mt')

    mt_final = mt.annotate_cols(
        additive_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].beta, 5, ordering = mt.additive[0].pval)),
        dominance_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].beta, 5, ordering = mt.additive[0].pval)),
        additive_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].pval, 5, ordering = mt.additive[0].pval)),
        dominance_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].pval, 5, ordering = mt.additive[0].pval)), 
        AC = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].AC, 5, ordering = mt.additive[0].pval)),
        n_complete_samples = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].n_complete_samples, 5, ordering = mt.additive[0].pval)) 
    )  

    mt_final = mt_final.annotate_cols(MAF = (mt_final.AC / (2 * mt_final.n_complete_samples)))

    mt_final = mt_final.annotate_cols(additive_hits_check = hl.map(
        lambda i: hl.sqrt(2 * mt_final.MAF[i] * (1.0 - mt_final.MAF[i])),
        [0,1,2,3,4]))
    mt_final = mt_final.annotate_cols(additive_hits_rescaled = mt_final.additive_hits / mt_final.additive_hits_check)

    cols = mt_final.cols()
    cols = cols.annotate(**{f'beta_A_{i}': cols.additive_hits_rescaled[i] for i in range(5)})
    cols = cols.annotate(**{f'beta_D_{i}': cols.dominance_hits[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_A_{i}': cols.additive_pval[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_D_{i}': cols.dominance_pval[i] for i in range(5)})

    cols.drop('additive_hits', 'dominance_hits', 'additive_pval', 'dominance_pval').export(f'{main_gs_folder}/{file_root_tsv}additive_maf_{MAF_query}.tsv.bgz')

# Next, the same but at dominance loci

# April 2021.
# Next, do the same but filter to dominance hits and consider the additive hits at these same locations.
# Fig. 3 - alternate version, for supplement.
mt = hl.read_matrix_table(f'{main_gs_folder}/{file_root}.mt')

v = hl.import_table(f'{main_gs_folder}/all_mafge0.05.hg19_multianno.csv', delimiter=',', quote = '"', missing='.', impute=True)
v = v.annotate(locus = hl.locus(contig = v.Chr, pos = v.Start), alleles = [v.Ref, v.Alt])
v = v.key_by(v.locus, v.alleles)

v = v.select(v.cytoBand)

mt = mt.annotate_rows(annovar = v[mt.row_key])
mt = mt.annotate_entries(MAF = (mt.additive.AC / (2 * mt.additive.n_complete_samples)))
mt = mt.annotate_entries(MAF = hl.case().when(mt.MAF > 0.5, 1 - mt.MAF)
                                                    .default(mt.MAF))
mt = mt.filter_entries(mt.additive.low_confidence_variant == False)
mt = mt.filter_rows(mt.locus.in_autosome())
mt = mt.checkpoint(f'{main_gs_folder}/group_by_cytoband/{file_root}_dominance_hits_by_cytoband_checkpoint.mt', overwrite=True)

for MAF_query in [0.01, 0.05]:
    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_hits_by_cytoband_checkpoint.mt')
    if str(MAF_query) == '0.01':
        mt = mt.filter_rows(mt.p_hwe > 1e-10)
    else:
        mt = mt.filter_rows(mt.p_hwe > 1e-6)
    mt.group_rows_by(mt.annovar.cytoBand).aggregate(
        additive = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.additive, 1, ordering = mt.dominance.pval)),
        dominance = hl.agg.filter(mt.MAF > MAF_query, hl.agg.take(mt.dominance, 1, ordering = mt.dominance.pval))).write(f'{main_gs_folder}/group_by_cytoband/{file_root}_dominance_hits_by_cytoband.mt', overwrite=True)

    mt = hl.read_matrix_table(f'{main_gs_folder}/group_by_cytoband/{file_root}_dominance_hits_by_cytoband.mt')

    mt_final = mt.annotate_cols(
        additive_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].beta, 5, ordering = mt.dominance[0].pval)),
        dominance_hits = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].beta, 5, ordering = mt.dominance[0].pval)),
        additive_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                      hl.agg.take(mt.additive[0].pval, 5, ordering = mt.dominance[0].pval)),
        dominance_pval = hl.agg.filter(hl.len(mt.additive) > 0,
                                       hl.agg.take(mt.dominance[0].pval, 5, ordering = mt.dominance[0].pval)), 
        AC = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].AC, 5, ordering = mt.dominance[0].pval)),
        n_complete_samples = hl.agg.filter(hl.len(mt.additive) > 0,
                            hl.agg.take(mt.additive[0].n_complete_samples, 5, ordering = mt.dominance[0].pval)) 
    )  

    mt_final = mt_final.annotate_cols(MAF = (mt_final.AC / (2 * mt_final.n_complete_samples)))

    mt_final = mt_final.annotate_cols(additive_hits_check = hl.map(
        lambda i: hl.sqrt(2 * mt_final.MAF[i] * (1.0 - mt_final.MAF[i])),
        [0,1,2,3,4]))
    mt_final = mt_final.annotate_cols(additive_hits_rescaled = mt_final.additive_hits/mt_final.additive_hits_check)

    cols = mt_final.cols()
    cols = cols.annotate(**{f'beta_A_{i}': cols.additive_hits_rescaled[i] for i in range(5)})
    cols = cols.annotate(**{f'beta_D_{i}': cols.dominance_hits[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_A_{i}': cols.additive_pval[i] for i in range(5)})
    cols = cols.annotate(**{f'pval_D_{i}': cols.dominance_pval[i] for i in range(5)})

    cols.drop('additive_hits', 'dominance_hits', 'additive_pval', 'dominance_pval').export(f'{main_gs_folder}/{file_root_tsv}dominance_no_X_maf_{MAF_query}.tsv.bgz')

# Finally, do the same as above, but enforce that both are genome-wide significant. If they're not, set to 0. This is done in the plotting script using R.

