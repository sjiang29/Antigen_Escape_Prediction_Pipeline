#!/bin/bash
# Generate DMS files (mutants + WT) and a cleaned A3M for Tranception

FASTA="chikv_e1_from_e1e2e3.fasta"
MSA_RAW="chikv_e1_from_e1e2e3.a3m"
OUT_MUTANTS="chikv_e1_from_e1e2e3_dms_mutants.csv"
OUT_WT="chikv_e1_from_e1e2e3_dms_wt.csv"
OUT_REPAIRED="chikv_e1_from_e1e2e3_ready.a3m"

python make_dms_and_fix_a3m.py \
  "$FASTA" \
  "$MSA_RAW" \
  "$OUT_MUTANTS" \
  "$OUT_WT" \
  "$OUT_REPAIRED"
