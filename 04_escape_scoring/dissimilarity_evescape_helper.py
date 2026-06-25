#!/usr/bin/env python3
"""
dissimilarity_evescape_helper.py
--------------------------------
EVEscape-style dissimilarity (p_d) from a mutants list.

Implements the paper's metric:
- Compute Δcharge and Δhydropathy between WT and mutant residue
  (charge in {-1,0,+1} at ~pH 7; hydropathy = Eisenberg–Weiss consensus scale).
- Standard-scale each difference across the protein's mutants and SUM the two standardized values.
- Temperature-scaled sigmoid to map to (0,1) => p_d.

Inputs
  --mutants_csv MUTANTS_CSV  : CSV with 'mutant' column (A123V)
  --out_csv OUT_CSV          : Output CSV with columns: mutant, p_d
  --tau TAU                  : Temperature for sigmoid (default 2.0)
  --start_idx START_IDX      : Optional inclusive WT FASTA start position
  --end_idx END_IDX          : Optional inclusive WT FASTA end position

Notes:
- If start_idx/end_idx are provided, mutants outside the requested WT sequence
  range are filtered out BEFORE standardization and p_d calculation.
"""
import argparse
import re
import numpy as np
import pandas as pd
from scipy.special import expit

HYDRO = {
    'A': 0.62, 'R': -2.53, 'N': -0.78, 'D': -0.90, 'C': 0.29, 'Q': -0.85, 'E': -0.74,
    'G': 0.48, 'H': -0.40, 'I': 1.38, 'L': 1.06, 'K': -1.50, 'M': 0.64, 'F': 1.19,
    'P': 0.12, 'S': -0.18, 'T': -0.05, 'W': 0.81, 'Y': 0.26, 'V': 1.08
}

CHARGE = {
    'D': -1, 'E': -1, 'K': +1, 'R': +1, 'H': +1,
    'A': 0, 'C': 0, 'F': 0, 'G': 0, 'I': 0, 'L': 0, 'M': 0, 'N': 0, 'Q': 0,
    'P': 0, 'S': 0, 'T': 0, 'V': 0, 'W': 0, 'Y': 0
}


def parse_mut(label):
    m = re.match(r'^([A-Z])(\d+)([A-Z])$', str(label))
    return (m.group(1), int(m.group(2)), m.group(3)) if m else (None, None, None)


def validate_range_from_df(start_idx, end_idx, parsed_positions):
    if start_idx is None:
        start_idx = int(parsed_positions.min())
    if end_idx is None:
        end_idx = int(parsed_positions.max())

    if start_idx < 1 or end_idx < 1:
        raise ValueError("start_idx and end_idx must be >= 1.")
    if start_idx > end_idx:
        raise ValueError(f"start_idx ({start_idx}) cannot be greater than end_idx ({end_idx}).")

    return start_idx, end_idx


def main():
    ap = argparse.ArgumentParser(description="EVEscape-style dissimilarity (p_d)")
    ap.add_argument("--mutants_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--tau", type=float, default=2.0)
    ap.add_argument("--start_idx", type=int, default=None,
                    help="Inclusive WT FASTA start position")
    ap.add_argument("--end_idx", type=int, default=None,
                    help="Inclusive WT FASTA end position")
    args = ap.parse_args()

    df = pd.read_csv(args.mutants_csv)
    if "mutant" not in df.columns:
        raise SystemExit("mutants_csv must contain a 'mutant' column.")

    parsed_rows = []
    for lab in df["mutant"].astype(str):
        w, pos, v = parse_mut(lab)
        if w is None or w not in HYDRO or v not in HYDRO:
            continue
        parsed_rows.append({"mutant": lab, "wt": w, "pos": pos, "mut": v})

    parsed = pd.DataFrame(parsed_rows)
    if parsed.empty:
        raise SystemExit("No valid mutants parsed.")

    start_idx, end_idx = validate_range_from_df(args.start_idx, args.end_idx, parsed["pos"])

    parsed = parsed[(parsed["pos"] >= start_idx) & (parsed["pos"] <= end_idx)].copy()
    if parsed.empty:
        raise SystemExit(
            f"No valid mutants remain after filtering to WT range {start_idx}-{end_idx}."
        )

    parsed["d_h"] = parsed.apply(lambda r: abs(HYDRO[r["wt"]] - HYDRO[r["mut"]]), axis=1)
    parsed["d_c"] = parsed.apply(lambda r: abs(CHARGE.get(r["wt"], 0) - CHARGE.get(r["mut"], 0)), axis=1)

    for col in ["d_h", "d_c"]:
        mu = float(parsed[col].mean())
        sd = float(parsed[col].std(ddof=0)) or 1.0
        parsed[col + "_z"] = (parsed[col] - mu) / sd

    parsed["dissim_combined"] = parsed["d_h_z"] + parsed["d_c_z"]

    mu = float(parsed["dissim_combined"].mean())
    sd = float(parsed["dissim_combined"].std(ddof=0)) or 1.0
    z = (parsed["dissim_combined"] - mu) / sd
    parsed["p_d"] = np.clip(expit(z / float(args.tau)), 1e-9, 1 - 1e-9)

    out = parsed[["mutant", "p_d"]].drop_duplicates("mutant")
    out.to_csv(args.out_csv, index=False)
    print(f"[dissimilarity] WT range {start_idx}-{end_idx}; wrote {len(out)} rows to {args.out_csv}")


if __name__ == "__main__":
    main()
