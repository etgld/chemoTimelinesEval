# Melanoma
## Strict
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/melanoma_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/melanoma.tsv \
       --output_dir subtask1_test_error_analysis/melanoma/strict/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/melanoma_test_patient_ids.txt \
       --strict
## Day
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/melanoma_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/melanoma.tsv \
       --output_dir subtask1_test_error_analysis/melanoma/day/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/melanoma_test_patient_ids.txt \
       --relaxed_to day
## Month
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/melanoma_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/melanoma.tsv \
       --output_dir subtask1_test_error_analysis/melanoma/month/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/melanoma_test_patient_ids.txt \
       --relaxed_to month
## Year
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/melanoma_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/melanoma.tsv \
       --output_dir subtask1_test_error_analysis/melanoma/year/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/melanoma_test_patient_ids.txt \
       --relaxed_to year
# Breast
## Strict
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/breast_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/breast.tsv \
       --output_dir subtask1_test_error_analysis/breast/strict/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/breast_test_patient_ids.txt \
       --strict
## Day
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/breast_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/breast.tsv \
       --output_dir subtask1_test_error_analysis/breast/day/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/breast_test_patient_ids.txt \
       --relaxed_to day
## Month
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/breast_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/breast.tsv \
       --output_dir subtask1_test_error_analysis/breast/month/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/breast_test_patient_ids.txt \
       --relaxed_to month
## Year
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/breast_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/breast.tsv \
       --output_dir subtask1_test_error_analysis/breast/year/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/breast_test_patient_ids.txt \
       --relaxed_to year
# Ovarian
## Strict
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/ovarian_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/ovarian.tsv \
       --output_dir subtask1_test_error_analysis/ovarian/strict/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/ovarian_test_patient_ids.txt \
       --strict
## Day
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/ovarian_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/ovarian.tsv \
       --output_dir subtask1_test_error_analysis/ovarian/day/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/ovarian_test_patient_ids.txt \
       --relaxed_to day
## Month
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/ovarian_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/ovarian.tsv \
       --output_dir subtask1_test_error_analysis/ovarian/month/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/ovarian_test_patient_ids.txt \
       --relaxed_to month
## Year
python summarize_evaluate_analyze.py --gold_path ~/HDD/DockerReadyTimelinesData/test/gold_timelines/ovarian_test_anafora_gold_timelines.json \
       --input_tsv ../deepphe_timelines_open/subtask1_test_timelines/ovarian.tsv \
       --output_dir subtask1_test_error_analysis/ovarian/year/ \
       --all_id_path ~/HDD/DockerReadyTimelinesData/test/all_patient_notes/patient_ids/ovarian_test_patient_ids.txt \
       --relaxed_to year
