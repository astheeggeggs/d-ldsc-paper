#!/bin/bash

gsutil hash -h gs://ukb-mega-gwas-results/round2/additive-tsvs/*tsv.bgz > md5sum_cloud.tsv

# Now cleanup in R to get the md5sums for comparison.
# I do this in the file additive_md5sum_cloud_cleanup.r