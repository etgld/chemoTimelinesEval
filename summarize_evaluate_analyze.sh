#!/bin/bash

# arg 1 is docker tsv
# arg 2 is cancer type/cohort
# arg 3 is split (train/dev/test)
# arg 4 gold data dir
# arg 5 eval mode
TMPDIR=".tmp_${2}_${3}_$(date +%F_%T)"
ERRDIR="error_analysis_${2}_${3}_$(date +%F_%T)"
mkdir $TMPDIR
mkdir $ERRDIR

python docker_output_to_timeline.py --input_tsv $1 --cancer_type $2 --output_dir $TMPDIR
cd $TMPDIR
SYSTEMTIMELINE="${2}_${3}_system_timelines.json"

case $5 in

    strict)
        python eval_timeline.py --pred_path $SYSTEMTIMELINE
        ;;
    day)
        python eval_timeline.py --pred_path $SYSTEMTIMELINE
        ;;
    month)
        python eval_timeline.py --pred_path $SYSTEMTIMELINE
        ;;
    year)
        python eval_timeline.py --pred_path $SYSTEMTIMELINE
        ;;
    *)
        echo -n "Incorrect mode, choose one of strict, day, month, year"
        ;;
    esac
