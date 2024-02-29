#!/bin/bash

# arg 1 is docker tsv
# arg 2 is cancer type/cohort
# arg 3 is split (train/dev/test)
# arg 4 gold data dir
# arg 5 eval mode
DOCKERTSV="$(readlink -m $1)"
CTYPE=$2
MODEARG=$5
TMPDIR="$(readlink -m ./.tmp_${2}_${3}_$(date +%a_%b_%d_%Y_%H_%M_%S_%p_%Z))"
ERRDIR="$(readlink -m ./error_analysis_${2}_${3}_$(date +%a_%b_%d_%Y_%H_%M_%S_%p_%Z))"
TASKDIR="$(readlink -m ${4}/subtask1/)"
# shouldn't need readlink for these since TASKDIR
# should already be absolute
ALLIDS="${TASKDIR}/All_Patient_IDs/${2}_${3}_patient_ids.txt"
GOLDTIMELINE="${TASKDIR}/Gold_Timelines_allPatients/${2}_${3}_all_patients_gold_timelines.json"
mkdir $TMPDIR
mkdir $ERRDIR

python docker_output_to_timeline.py --input_tsv $DOCKERTSV --cancer_type $CTYPE --output_dir $TMPDIR
cd $TMPDIR
# SYSTEMTIMELINE="${2}_${3}_system_timelines.json"
SYSTEMTIMELINE="${2}_timelines.json"

case $MODEARG in
    strict)
        EVALMODE="--strict"
        ;;
    day)
        EVALMODE="--relaxed_to day"
        ;;
    month)
        EVALMODE="--relaxed_to month"
        ;;
    year)
        EVALMODE="--relaxed_to year"
        ;;
    *)
        echo -n "Incorrect mode, choose one of strict, day, month, year"
        exit 1
        ;;
esac

python ../eval_timeline.py --debug --pred_path $SYSTEMTIMELINE --gold_path $GOLDTIMELINE --all_id_path $ALLIDS $EVALMODE
python ../docker_output_error_analysis.py --docker_tsv $DOCKERTSV --error_json ./patient_level_debug.json --output_dir $ERRDIR $EVALMODE
