#!/bin/bash

# Biomarkers
for pipeline in {3..6}
do
   hailctl dataproc submit dp export_dominance_results_biomarkers_no_HWE.py $pipeline
   hailctl dataproc submit dp export_dominance_results_biomarkers_no_HWE_GP.py $pipeline
done