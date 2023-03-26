#\!/bin/bash

source /broad/software/scripts/useuse

export GOOGLE_APPLICATION_CREDENTIALS=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/inputs/ukbb-round2-a700709c0e9a.json
use BEDTools

# Ensure redhat 7 compute nodes.
# Get the environment that was used to generate these data.
# You must clone and run Masa's script to run the finemapping code 
# https://github.com/mkanai/xfinemap

# The following needs to be on the cluster - this is what I need to loop over to create the 
# finemapping inputs

# pip install pybedtools --user
# pip install python-dateutil --user
# pip install bedtools --user
# pip install htslib --user
# pip install google-auth --user
# pip install dsub --user

pheno=$1
# phenotype name (this is looped over in the submit_all_make_finemap_inputs_hail_input.r script).
full_phenotype_file=$2
# File from which you can extract this phenotype in $1.
# This is one of:
# 1. /psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping/ukb30163.restricted_and_phesant_formatted.biomarkers.tsv.gz
# 2. /psych/genetics_data/projects/ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz
# 3. /psych/genetics_data/projects/ukb31063/ukb31063.ICD10_phenotypes.both_sexes.tsv.bgz
# 4. /psych/genetics_data/projects/ukb31063/ukb31063.FINNGEN_phenotypes.both_sexes.tsv.bgz

maf_locus=$7
# I ran 0.05.
# maf_locus=0.05

maf_finemap=$8
# maf_finemap=0
# This is the MAF cutoff for the phenotype. After merging the additive and dominance summary 
# statistics files, the variants are filtered > ${maf_locus}

p_hwe_finemap=$9
# This is the p_HWE cutoff. After merging the additive and dominance summary statistics files,
# the variants are filtered > ${p_hwe_finemap}
# p_hwe_finemap=1e-10

incl_file=$3 # Note that the incl file is on the cluster and is transferred to the cloud later in the script.
# The location of the .incl file
# This has format:
# /psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping/inputs/incl_files/{phesant,icd10,finngen,biomarkers}/${pheno}.incl

dom_file=$4
# The location of the dominance summary statistics file.
# This has format:
# /psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/both_sexes/biomarkers/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz for biomarkers
# and /psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/both_sexes/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz otherwise.

add_file=$5
# The location of the additive summary statistics file.
# This has format:
# /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/biomarkers/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz for biomarkers
# and /stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz otherwise.

out=$6

# The location of the output of this initial script: 
# This has format:
# /psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping/{phesant,icd10,finngen,biomarkers}
# 
# This folder will contain the following collection of subfolders.
# A. ${out}/inputs/both_sexes_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
#   1. The merged finemapping input file: ${pheno}_finemapping_input.gz
#   2. The .json file for additive finemapping: ${pheno}_dom.json.
#   3. The .json file for dominance finemapping: ${pheno}_add.json.
# B. ${out}/inputs/z_files_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
#   1. ${pheno}_finemapping_input.chr{contig}.{region_start}-{region_end}.
#   2. ${pheno}_finemapping_input.finemap.dsub.sh, ${pheno}_finemapping_input.finemap.tasks.txt
#   3. ${pheno}_finemapping_input.ldstore.dsub.sh, ${pheno}_finemapping_input.ldstore.tasks.txt
#   4. ${pheno}_finemapping_input.susie.dsub.sh, ${pheno}_finemapping_input.susie.tasks.txt
# C. ${out}/finemapping_inputs/both_sexes_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
#   1. ${pheno}.bed: tab delimited file with open intervals of significant dominance loci.
#   2. ${pheno}_dom.log: log file of the dominance call to ${script} (/home/unix/dpalmer/Repositories/xfinemap/make_finemap_inputs.py).
#   3. ${pheno}_add.log: log file of the additive call to ${script} (/home/unix/dpalmer/Repositories/xfinemap/make_finemap_inputs.py).
#   4. ${pheno}.lead_snps.txt: The collection of lead snps that are greedily added to expand the locus.

# Some examples to use for testing.
# cat example
# pheno=1747_2
# full_phenotype_file=/psych/genetics_data/projects/ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz
# pheno=1150_3
# maf_locus=0.05
# maf_finemap=0
# p_hwe_finemap=1e-10
# dom_file=/psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/both_sexes/${pheno}.dominance.gwas.imputed_v3.both_sexes.tsv.bgz
# add_file=/stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz
# out=/psych/genetics_data/dpalmer/UKbb/testing/
# incl_file=/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping/inputs/incl_files/phesant/${pheno}.incl

# cts example
# full_phenotype_file=/psych/genetics_data/projects/ukb31063/ukb31063.PHESANT_January_2019.both_sexes.tsv.bgz
# pheno=50_irnt
# maf_locus=0.05
# maf_finemap=0
# p_hwe_finemap=1e-10
# dom_file=/psych/genetics_data/dpalmer/UKbb/sumstats-dominance-export/both_sexes/${pheno}.dominance.gwas.imputed_v3.both_sexes.tsv.bgz
# add_file=/stanley/genetics/analysis/ukbb_sumstats/sumstats-additive-export/both_sexes/${pheno}.gwas.imputed_v3.both_sexes.tsv.bgz
# out=/psych/genetics_data/dpalmer/UKbb/testing/
# incl_file=/psych/genetics_data/dpalmer/UKbb/dominance_fine_mapping/inputs/incl_files/phesant/${pheno}.incl

outpath=${out}/inputs/both_sexes_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
mkdir -p $outpath

finemapping_outpath=$out/finemapping_inputs/both_sexes_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
finemapping_outpath_dom=${finemapping_outpath}/dominance
finemapping_outpath_add=${finemapping_outpath}/additive

mkdir -p $finemapping_outpath_dom
mkdir -p $finemapping_outpath_add

script=/home/unix/dpalmer/Repositories/xfinemap/make_finemap_inputs.py
json_script=/home/unix/dpalmer/Repositories/ldscgxe/dominance_fine_mapping/create_json_files.r

# Note, creates the z files in the local directory and then pushes them to the cloud.
zout=${out}/inputs/z_and_task_files_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}
zout_dom=${zout}/dominance
zout_add=${zout}/additive

mkdir -p $zout_dom
mkdir -p $zout_add

gsout_add=gs://ukbb-dominance/dominance_fine_mapping_rerun/additive_output_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}/
gsout_dom=gs://ukbb-dominance/dominance_fine_mapping_rerun/dominance_output_maf_${maf_locus}_p_hwe_${p_hwe_finemap}_maf_finemap_${maf_finemap}/

incl_file_cloud="gs://ukbb-dominance/dominance_fine_mapping_rerun/output/incl_files/${pheno}.incl"
gsutil cp "${incl_file}" "${incl_file_cloud}"
input_samples="gs://ukbb-dominance/dominance_fine_mapping/ukb31063.sample"
jsonfile_loc=${outpath}/${pheno}

Rscript $json_script $pheno $gsout_add $gsout_dom $incl_file $incl_file_cloud $n_samples $input_samples $jsonfile_loc

jsonfile_dom_loc="${jsonfile_loc}_dom.json"
jsonfile_add_loc="${jsonfile_loc}_add.json"

variants_mapping_file=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/variants.tsv.bgz

add_template_head=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/additive_header.txt
add_template_head_ord=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/additive_header_ord.txt
dom_template_head=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/dominance_header.txt

# There are different columns in cts vs cat variables - important!
add_template_head_cts=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/additive_header_cts.txt
dom_template_head_cts=/stanley/genetics_storage/analysis/ukbb_dominance/dominance_fine_mapping/dominance_header.txt # It's the same for the dominance results.

finemapping_input_head=/stanley/genetics/analysis/ukbb_dominance/results/gwas/dominantSig/output/merged/genotyped/header_withMAF_for_finemapping_duncan.txt
finemapping_input_head_both=/stanley/genetics/analysis/ukbb_dominance/results/gwas/dominantSig/output/merged/genotyped/header_withMAF_for_finemapping_add_and_dom_duncan.txt

if (diff <(head -1 $add_template_head) <(zcat $add_file | head -1)) || (diff <(head -1 $add_template_head_ord) <(zcat $add_file | head -1)); then 
    echo "additive header is equal";
    if diff <(head -1 $dom_template_head) <(zcat $dom_file | head -1); then 
        
        echo "dominance header is equal";

        if [ -f "$outpath/${pheno}_finemapping_input.gz" ]; then
            echo "$outpath/${pheno}_finemapping_input.gz exists, don't need to create merged summary statistics file"
        else 
            echo "$outpath/${pheno}_finemapping_input.gz does not exist, creating merged summary statistics file."

            (cat ${finemapping_input_head_both} \
                <(
                    (
                        (LANG=en_EN join -t $'\t' -j 1 -o 1.1,2.2,2.3,2.9,2.10,2.11,2.12,1.4,1.5,1.6,1.7,2.9,2.10,2.11,2.12 <(zcat "${dom_file}" | tail -n +2 | LANG=en_EN sort -k 1,1) <(zcat "${add_file}" | tail -n +2 | LANG=en_EN sort -k 1,1) | \
                            awk -F '[:\t]' '{print $1":"$2":"$3":"$4"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}') | \
                        LANG=en_EN join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.6,2.7,2.16 - <(zcat "${variants_mapping_file}" | tail -n +2 | LANG=en_EN sort -k 1,1)
                    ) | awk -F '\t' "\$7>${maf_finemap} && \$18>${p_hwe_finemap}" | \
                awk -F '\t' 'BEGIN { OFS="\t" } { if ($5==$6) $18=$7; else $18=1-$7; print $0; }'
                )
            ) | gzip > $outpath/${pheno}_finemapping_input.gz

        fi

        # What does this do? 
        # It joins based on the first column - if either of the columns are not ordered it will return an error, 
        # then we remove the first row, split on tab or : to then get the chromosome, position, ref and alt information,
        # the minor allele and the minor allele frequency, and the dominance p-value, se and beta information.
        # We then switch MAF to the reference allele frequency by checking whether the minor allele is equal to the 
        # reference allele using awk, finally we add the header. This then gets saved in a separate folder.

        # create (if required) and move to a new directory which will 
        # contain the task files for the dominance filemapping.

        cd ${zout_dom}

        python $script \
        --sumstats $outpath/${pheno}_finemapping_input.gz \
        --json ${jsonfile_dom_loc} \
        --out ${finemapping_outpath_dom}/${pheno}_dom \
        --exclude-MHC \
        --MHC-start 25000000\
        --dominant \
        --maf-threshold ${maf_locus}

        # create (if required) and move to a new directory which will 
        # contain the task files for the additive filemapping.
        
        cd ${zout_add}

        python $script \
        --sumstats ${outpath}/${pheno}_finemapping_input.gz \
        --json ${jsonfile_add_loc} \
        --out ${finemapping_outpath_add}/${pheno}_add \
        --exclude-MHC \
        --MHC-start 25000000\
        --maf-threshold ${maf_locus} \
        --bed ${finemapping_outpath_dom}/${pheno}_dom.bed \
        --no-upload

    else
        echo "dominance header is not equal, something's wrong.";
    fi
else 
    echo "additive header is not equal, hopefully this is a continuous variable"; 
    if diff <(head -1 $add_template_head_cts) <(zcat $add_file | head -1); then 
        echo "additive cts header is equal...looking good";
        if diff <(head -1 $dom_template_head_cts) <(zcat $dom_file | head -1); then 
    
            echo "dominance cts header is equal, great!";

            if [ -f "$outpath/${pheno}_finemapping_input.gz" ]; then
                echo "$outpath/${pheno}_finemapping_input.gz exists, don't need to create merged summary statistics file"
            else 
                echo "$outpath/${pheno}_finemapping_input.gz does not exist, creating merged summary statistics file."

                (cat ${finemapping_input_head_both} \
                    <(
                        (
                            (LANG=en_EN join -t $'\t' -j 1 -o 1.1,2.2,2.3,2.8,2.9,2.10,2.11,1.4,1.5,1.6,1.7,2.8,2.9,2.10,2.11 <(zcat "${dom_file}" | tail -n +2 | LANG=en_EN sort -k 1,1) <(zcat "${add_file}" | tail -n +2 | LANG=en_EN sort -k 1,1) | \
                                awk -F '[:\t]' '{print $1":"$2":"$3":"$4"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18}') | \
                            LANG=en_EN join -t $'\t' -j 1 -o 1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.6,2.7,2.16 - <(zcat "${variants_mapping_file}" | tail -n +2 | LANG=en_EN sort -k 1,1)
                        ) | awk -F '\t' "\$7>${maf_finemap} && \$18>${p_hwe_finemap}" | \
                    awk -F '\t' 'BEGIN { OFS="\t" } { if ($5==$6) $18=$7; else $18=1-$7; print $0; }'
                    )
                ) | gzip > ${outpath}/${pheno}_finemapping_input.gz
            fi

            # What does this do? 
            # It joins based on the first column - if either of the columns are not ordered it will return an error, 
            # then we remove the first row, split on tab or : to then get the chromosome, position, ref and alt information,
            # the minor allele and the minor allele frequency, and the dominance p-value, se and beta information.
            # We then switch MAF to the reference allele frequency by checking whether the minor allele is equal to the 
            # reference allele using awk, finally we add the header. This then gets saved in a separate folder.

            # create (if required) and move to a new directory which will 
            # contain the task files for the dominance filemapping.

            cd ${zout_dom}

            python $script \
            --sumstats ${outpath}/${pheno}_finemapping_input.gz \
            --json ${jsonfile_dom_loc} \
            --out ${finemapping_outpath_dom}/${pheno}_dom \
            --exclude-MHC \
            --MHC-start 25000000\
            --dominant \
            --maf-threshold ${maf_locus}

            # create (if required) and move to a new directory which will 
            # contain the task files for the additive filemapping.

            cd ${zout_add}

            python $script \
            --sumstats $outpath/${pheno}_finemapping_input.gz \
            --json ${jsonfile_add_loc} \
            --out ${finemapping_outpath_add}/${pheno}_add \
            --exclude-MHC \
            --MHC-start 25000000 \
            --maf-threshold ${maf_locus} \
            --bed ${finemapping_outpath_dom}/${pheno}_dom.bed \
            --no-upload

        fi
    fi
fi
