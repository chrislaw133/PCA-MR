#!/bin/bash

python run_QTL_analysis.py \
  --plink /path/to/plinkbinary \
  --annotation_file ${FEATURE_FILE} \
  --phenotype_file /path/to/cumulative-expression \
  --output_directory /path/to/outdir \
  --window 250000 \
  --covariates_file /path/to/covariates \
  --randomeff_files /path/to/kinship-matrix \
  --sample_mapping_file /path/to/sampling-mapping-file \
  --cis \
  --number_of_permutations 0 \
  --gaussianize_method ranknorm \
  --regress_covariates \
  --write_zscores
