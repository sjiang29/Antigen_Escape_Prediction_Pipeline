python evescape_pipeline.py \
  --wt_fasta EEEV_p62-E1_monomer_A.fasta \
  --wt_pdb EEEV_p62-E1_trimer.pdb \
  --tr_csv eeev_e1_from_e1e2e3_dms_mutants_0.8.csv \
  --out_dir pipeline_results_eeev_e1_from_e1e2e3 \
  --tau_f 1.0 --tau_a 1.0 --tau_d 2.0 \
  --chain A \
  --start_idx 447 \
  --end_idx 855
