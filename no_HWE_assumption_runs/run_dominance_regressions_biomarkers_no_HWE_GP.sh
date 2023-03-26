#!/bin/bash

# Biomarkers
for pipeline in {5..6}
do
   hailctl dataproc submit dp run_dominance_regressions_biomarkers_no_HWE_GP.py $pipeline
done