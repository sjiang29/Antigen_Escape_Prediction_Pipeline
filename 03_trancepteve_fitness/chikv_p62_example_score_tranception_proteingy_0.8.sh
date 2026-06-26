#!/bin/bash
# Run Tranception scoring for E2 mutants

# --- config ---
LEN=431  # replace with your WT length
FASTA="p62_chikv.fasta"
MSA="p62_chikv_ready.a3m"
DMS="p62_chikv_dms_mutants.csv"
OUTDIR="results_p62_chikv"
CHECKPOINT="Tranception_Large"

# --- run ---
python score_tranception_proteingym.py \
  --checkpoint $CHECKPOINT \
  --model_framework pytorch \
  --target_seq "$(awk 'NR>1{gsub(/[ \t\r\n]/,""); printf "%s",$0}' $FASTA)" \
  --MSA_folder . \
  --MSA_filename $MSA \
  --output_scores_folder $OUTDIR \
  --inference_time_retrieval \
  --retrieval_inference_weight 0.8 \
  --batch_size_inference 16 \
  --num_workers 4 \
  --DMS_data_folder . \
  --DMS_file_name $DMS \
  --MSA_start 1 \
  --MSA_end $LEN

