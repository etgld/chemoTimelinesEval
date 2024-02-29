#!/bin/bash

# arg 1 is docker tsv
# arg 2 is cancer type/cohort
# arg 3 is split (train/dev/test)

TMPDIR=".tmp_${2}_${3}_$(date +%F_%T)"
mkdir $TMPDIR

python docker_output_to_timeline.py --input_tsv $1 --cancer_type $2 --output_dir $TMPDIR
