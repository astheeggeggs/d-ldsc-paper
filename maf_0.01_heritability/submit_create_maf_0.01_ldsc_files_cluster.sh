#!/bin/bash

# Run it

for b in {1..31}; do
	qsub -v batch="$b" create_maf_0.01_ldsc_files_cluster_template.sh
done
