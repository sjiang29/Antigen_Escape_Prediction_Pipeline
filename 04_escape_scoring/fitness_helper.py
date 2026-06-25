#!/usr/bin/env python3
"""
fitness_helper.py
------------------
Build the EVEscape fitness component (p_f) from TranceptEVE output.

Inputs
  --wt_fasta WT_FASTA         : FASTA file with a single wild-type sequence
  --tr_csv TR_CSV             : TranceptEVE table with columns:
                                mutated_sequence, avg_score
                                (or both avg_score_L_to_R & avg_score_R_to_L)
  --out_csv OUT_CSV           : Output CSV path with columns:
                                mutant, avg_score, dllr_tr, z_f, p_f
  --tau TAU                   : Temperature for sigmoid after z-score (default: 2.0)
  --sep SEP                   : Optional CSV separator override (auto-detect by default)
  --keep_multi                : Keep non-single mutants in the output (default: drop)
  --start_idx START_IDX       : Optional inclusive WT FASTA start position
  --end_idx END_IDX           : Optional inclusive WT FASTA end position

Behavior
  1) If mutated_sequence length == WT FASTA length:
       - compare against full WT sequence
       - mutant labels use native WT numbering
  2) Else, if mutated_sequence length == (end_idx - start_idx + 1):
       - compare against WT subsequence wt_seq[start_idx-1:end_idx]
       - mutant labels are automatically shifted into global WT numbering
         e.g. local A1V with start_idx=451 becomes A451V
  3) Else:
       - error out with a clear message

Notes
  - Uses raw avg_score (per-residue average log-likelihood) to compute ΔLLR vs WT/reference.
  - If WT/reference sequence is not present in the CSV, falls back to mean-centering.
  - Single-mutant labels look like A123V.
"""

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit


def read_fasta_one_seq(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                continue
            seq.append(line)
    if not seq:
        raise ValueError(f"No sequence found in FASTA: {path}")
    return "".join(seq).strip()


def validate_range(start_idx, end_idx, seq_len):
    if start_idx is None:
        start_idx = 1
    if end_idx is None:
        end_idx = seq_len

    if start_idx < 1 or end_idx < 1:
        raise ValueError("start_idx and end_idx must be >= 1.")
    if start_idx > end_idx:
        raise ValueError(f"start_idx ({start_idx}) cannot be greater than end_idx ({end_idx}).")
    if end_idx > seq_len:
        raise ValueError(f"end_idx ({end_idx}) exceeds WT FASTA length ({seq_len}).")

    return start_idx, end_idx


def infer_mut_label(mut_seq: str, ref_seq: str, position_offset: int = 0):
    """
    Compare mutated sequence against reference sequence and infer a single-mutant label.
    Returns e.g. A123V, where numbering is shifted by position_offset if needed.
    Returns None if:
      - lengths differ
      - sequence is WT
      - sequence has >1 mutation
    """
    if len(mut_seq) != len(ref_seq):
        return None

    diffs = [
        (i, w, m)
        for i, (w, m) in enumerate(zip(ref_seq, mut_seq), start=1)
        if w != m
    ]

    if len(diffs) != 1:
        return None

    i, w, m = diffs[0]
    global_i = i + position_offset
    return f"{w}{global_i}{m}"


def extract_mut_pos(mutant_label: str):
    m = re.match(r"^[A-Z](\d+)[A-Z]$", str(mutant_label))
    return int(m.group(1)) if m else None


def get_reference_sequence_and_offset(df, wt_seq, start_idx, end_idx):
    """
    Decide whether mutated_sequence is full-length or corresponds to the requested region.

    Returns
      ref_seq, position_offset, mode, dominant_len
    """
    seq_lengths = df["mutated_sequence"].astype(str).str.len()
    if seq_lengths.empty:
        sys.exit("[fitness][error] No mutated_sequence values found.")

    length_counts = seq_lengths.value_counts()
    dominant_len = int(length_counts.index[0])

    print(f"[fitness] Detected most common mutated_sequence length: {dominant_len}")

    wt_len = len(wt_seq)
    region_len = end_idx - start_idx + 1

    if dominant_len == wt_len:
        ref_seq = wt_seq
        position_offset = 0
        mode = "full_length"
        print("[fitness] mutated_sequence matches full WT FASTA length.")
        print("[fitness] Using full WT sequence as reference.")
        return ref_seq, position_offset, mode, dominant_len

    if dominant_len == region_len:
        ref_seq = wt_seq[start_idx - 1:end_idx]
        position_offset = start_idx - 1
        mode = "subsequence"
        print("[fitness] mutated_sequence matches the requested WT subsequence length.")
        print(f"[fitness] Using WT subsequence positions {start_idx}-{end_idx} as reference.")
        print(f"[fitness] Auto-inferred position offset: {position_offset}")
        return ref_seq, position_offset, mode, dominant_len

    sys.exit(
        "[fitness][error] mutated_sequence length "
        f"({dominant_len}) does not match either:\n"
        f"  - WT FASTA length: {wt_len}\n"
        f"  - requested region length: {region_len} (from {start_idx}-{end_idx})\n"
        "Please check whether wt_fasta and tr_csv correspond to the same sequence region."
    )


def main():
    ap = argparse.ArgumentParser(description="Build p_f from TranceptEVE raw avg_score")
    ap.add_argument("--wt_fasta", required=True)
    ap.add_argument("--tr_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--tau", type=float, default=2.0)
    ap.add_argument("--sep", default=None)
    ap.add_argument("--keep_multi", action="store_true")
    ap.add_argument("--start_idx", type=int, default=None,
                    help="Inclusive WT FASTA start position")
    ap.add_argument("--end_idx", type=int, default=None,
                    help="Inclusive WT FASTA end position")
    args = ap.parse_args()

    wt_seq = read_fasta_one_seq(args.wt_fasta)
    wt_len = len(wt_seq)
    start_idx, end_idx = validate_range(args.start_idx, args.end_idx, wt_len)

    print(f"[fitness] WT FASTA length: {wt_len} aa")
    print(f"[fitness] Requested WT range: {start_idx}-{end_idx}")

    # Read CSV
    df = pd.read_csv(args.tr_csv, sep=args.sep if args.sep else None, engine="python")

    if "mutated_sequence" not in df.columns:
        sys.exit("[fitness][error] Need 'mutated_sequence' column.")

    # Ensure avg_score exists
    if "avg_score" not in df.columns:
        if {"avg_score_L_to_R", "avg_score_R_to_L"}.issubset(df.columns):
            df["avg_score"] = df[["avg_score_L_to_R", "avg_score_R_to_L"]].mean(axis=1)
            print("[fitness][warn] 'avg_score' missing; averaged L->R and R->L.")
        else:
            sys.exit(
                "[fitness][error] Need 'avg_score' or both "
                "'avg_score_L_to_R' and 'avg_score_R_to_L'."
            )

    # Decide reference sequence + numbering offset automatically
    ref_seq, position_offset, mode, dominant_len = get_reference_sequence_and_offset(
        df, wt_seq, start_idx, end_idx
    )

    # Keep only rows whose sequence length matches the chosen reference length
    before_len_filter = len(df)
    df = df[df["mutated_sequence"].astype(str).str.len() == len(ref_seq)].copy()
    after_len_filter = len(df)
    if after_len_filter < before_len_filter:
        print(
            f"[fitness][warn] Dropped {before_len_filter - after_len_filter} rows "
            f"with mutated_sequence length != {len(ref_seq)}."
        )

    # Label mutants
    df["mutant"] = df["mutated_sequence"].astype(str).apply(
        lambda s: infer_mut_label(s, ref_seq, position_offset)
    )

    single_mask = df["mutant"].notna()
    dropped = (~single_mask).sum()
    if dropped and not args.keep_multi:
        print(f"[fitness][warn] Dropping {dropped} non-single mutants. Use --keep_multi to keep.")
        df = df[single_mask].copy()

    if df.empty:
        sys.exit("[fitness][error] No valid single mutants remain after mutant parsing.")

    # Filter to requested global WT range
    df["mut_pos"] = df["mutant"].apply(extract_mut_pos)
    before_range_filter = len(df)
    df = df[
        df["mut_pos"].notna() &
        (df["mut_pos"] >= start_idx) &
        (df["mut_pos"] <= end_idx)
    ].copy()
    print(f"[fitness] Kept {len(df)} / {before_range_filter} mutants after WT range filtering.")

    if df.empty:
        sys.exit(f"[fitness][error] No valid single mutants remain in WT range {start_idx}-{end_idx}.")

    # Find WT/reference score if present in same-length subset
    ref_rows = df.loc[df["mutated_sequence"] == ref_seq]
    if len(ref_rows) > 0:
        ref_avg = float(ref_rows["avg_score"].mean())
        if mode == "full_length":
            print(f"[fitness] Found WT in CSV. WT avg_score={ref_avg:.6f}")
        else:
            print(f"[fitness] Found reference subsequence in CSV. Reference avg_score={ref_avg:.6f}")
    else:
        ref_avg = float(df["avg_score"].mean())
        if mode == "full_length":
            print("[fitness][warn] WT not found in filtered CSV. Using mean(avg_score) for centering.")
        else:
            print("[fitness][warn] Reference subsequence not found in filtered CSV. Using mean(avg_score) for centering.")

    # ΔLLR per residue
    df["dllr_tr"] = df["avg_score"] - ref_avg

    # Standardize and map to (0,1)
    mu = float(df["dllr_tr"].mean())
    sd = float(df["dllr_tr"].std(ddof=0)) or 1.0
    df["z_f"] = (df["dllr_tr"] - mu) / sd
    df["p_f"] = np.clip(expit(df["z_f"] / float(args.tau)), 1e-9, 1 - 1e-9)

    # Keep tidy columns
    keep_cols = ["mutant", "mutated_sequence", "avg_score", "dllr_tr", "z_f", "p_f"]
    for extra in ["avg_score_L_to_R", "avg_score_R_to_L"]:
        if extra in df.columns:
            keep_cols.append(extra)

    out = (
        df[keep_cols]
        .dropna(subset=["mutant"])
        .groupby("mutant", as_index=False)
        .mean(numeric_only=True)
    )

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)

    mutants_path = str(Path(args.out_csv).with_suffix(".mutants.csv"))
    out[["mutant"]].to_csv(mutants_path, index=False)

    print(f"[fitness][done] Wrote {len(out)} single mutants to: {args.out_csv}")
    print(f"[fitness][done] Mutant list: {mutants_path}")


if __name__ == "__main__":
    main()
