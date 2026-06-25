#!/usr/bin/env python3
"""
accessibility_evescape_helper.py
--------------------------------
EVEscape-style accessibility (p_a) from a WT PDB + WT FASTA + mutants list.

Implements the paper's metric:
- Accessibility per residue = NEGATIVE weighted contact number (WCN)
- WCN_i = sum_{j != i} 1 / d_ij^2,
  where d_ij is the distance between side-chain geometric centers of residues i and j.
  (If a residue has no side chain atoms, fall back to Cα.)
- Aggregate across the whole biological assembly provided (all chains in the PDB model).
- After computing WCN, take NEGATIVE WCN as the accessibility signal,
  then per-protein z-score and temperature-scaled sigmoid to get p_a in (0,1).

Inputs
  --wt_fasta WT_FASTA            : FASTA with the WT sequence (1 sequence)
  --wt_pdb WT_PDB                : PDB of the WT protein (extracellular conformation)
  --mutants_csv MUTANTS_CSV      : CSV with a 'mutant' column (e.g., A123V)
  --out_csv OUT_CSV              : Output CSV with columns: mutant, p_a
  --tau TAU                      : Temperature for sigmoid (default 2.0)
  --chain CHAIN                  : Optional; if provided, restrict to a single chain ID
  --allow_gaps                   : If set, linearly impute missing positions after mapping
  --start_idx START_IDX          : Optional inclusive WT FASTA start position
  --end_idx END_IDX              : Optional inclusive WT FASTA end position

Dependencies: biopython, scipy
Notes:
- We align the PDB-derived sequence to the WT FASTA to map FASTA indices to PDB residues.
- If your PDB contains multiple chains of the same protomer (e.g., a trimer), we include ALL of
  them in WCN (as in the paper: “full Spike trimer”). If you want only one chain, pass --chain.
- start_idx/end_idx only restrict which WT FASTA positions are returned in the output.
  Structural WCN is still computed using the full selected structure context.
"""
import argparse
import re
import numpy as np
import pandas as pd
from Bio import PDB
from Bio.PDB.Polypeptide import PPBuilder
from Bio import pairwise2
from scipy.special import expit

BB_NAMES = {"N", "CA", "C", "O"}  # backbone atom names


def read_fasta_one_seq(path):
    seq = []
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith(">"):
                continue
            seq.append(line)
    if not seq:
        raise ValueError("No sequence found in FASTA.")
    return "".join(seq)


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


def pdb_chain_sequences(structure, chain_id=None):
    """Return list of (chain_id, sequence, residues) where residues is a list of PDB Residue objects
    matching the sequence (1 letter per standard AA)."""
    ppb = PPBuilder()
    results = []
    for model in structure:
        for ch in model:
            if chain_id and ch.id != chain_id:
                continue
            seq = ""
            res_list = []
            for pp in ppb.build_peptides(ch):
                seq += str(pp.get_sequence())
                res_list += pp
            if seq:
                results.append((ch.id, seq, res_list))
        break  # use first model
    if chain_id and not results:
        raise RuntimeError(f"Chain {chain_id} not found.")
    return results


def sidechain_center(res):
    """Geometric center of side-chain atoms; fallback to CA if no side-chain."""
    atoms = [a for a in res.get_atoms() if a.get_name() not in BB_NAMES]
    if not atoms:
        ca = res["CA"] if "CA" in res else None
        if ca is None:
            return None
        return np.array(ca.get_coord(), dtype=float)
    coords = np.array([a.get_coord() for a in atoms], dtype=float)
    return coords.mean(axis=0)


def compute_wcn_allchains(structure, chain_id=None):
    """Compute WCN per residue index across the provided chains (or a selected chain).
    Returns DataFrame with columns: chain, pdb_index (0..n-1 within that chain's polypeptide),
    resid (PDB number), wcn, and global_index."""
    chains = pdb_chain_sequences(structure, chain_id)

    centers = []
    meta = []  # (chain, pdb_index, resid, residue_object)
    for (cid, seq, residues) in chains:
        for i, res in enumerate(residues):
            cen = sidechain_center(res)
            if cen is None:
                cen = np.array([np.nan, np.nan, np.nan])
            centers.append(cen)
            meta.append((cid, i, res.id[1], res))

    centers = np.array(centers, dtype=float)
    n = len(centers)
    wcn_vals = np.zeros(n, dtype=float)

    for i in range(n):
        ci = centers[i]
        if np.any(np.isnan(ci)):
            wcn_vals[i] = np.nan
            continue
        diffs = centers - ci
        d2 = np.sum(diffs * diffs, axis=1)
        d2[i] = np.nan
        inv = 1.0 / d2
        inv[np.isinf(inv)] = np.nan
        wcn_vals[i] = np.nansum(inv)

    rows = []
    for idx, ((cid, pdb_index, resid, _), wcn) in enumerate(zip(meta, wcn_vals)):
        rows.append({
            "global_index": idx,
            "chain": cid,
            "pdb_index": pdb_index,
            "resid": resid,
            "wcn": wcn
        })
    return pd.DataFrame(rows)


def align_fasta_to_pdbseq(wt_fasta_seq, pdb_seq_concat):
    """Global alignment to map FASTA positions (1-based) -> PDB concat positions (1-based)."""
    alns = pairwise2.align.globalxx(wt_fasta_seq, pdb_seq_concat, one_alignment_only=True)
    if not alns:
        raise RuntimeError("Alignment failed between FASTA and PDB sequence.")
    a, b, *_ = alns[0]

    fasta_to_pdb_concat = {}
    fi = pi = 0
    pos_f = pos_p = 0
    while fi < len(a) and pi < len(b):
        if a[fi] != "-":
            pos_f += 1
        if b[pi] != "-":
            pos_p += 1
        if a[fi] != "-" and b[pi] != "-":
            fasta_to_pdb_concat[pos_f] = pos_p
        fi += 1
        pi += 1
    return fasta_to_pdb_concat


def concat_chain_sequences(structure, chain_id=None):
    chains = pdb_chain_sequences(structure, chain_id)
    concat = ""
    concat_to_global = {}

    df_all = compute_wcn_allchains(structure, chain_id)

    for (cid, seq, residues) in chains:
        for j, res in enumerate(residues):
            row = df_all[(df_all["chain"] == cid) & (df_all["pdb_index"] == j)].iloc[0]
            concat += seq[j]
            idx_concat = len(concat)  # 1-based
            concat_to_global[idx_concat] = int(row["global_index"])

    return concat, concat_to_global, df_all


def parse_pos(label):
    m = re.match(r"^[A-Z]([0-9]+)[A-Z]$", str(label))
    return int(m.group(1)) if m else None


def main():
    ap = argparse.ArgumentParser(description="EVEscape-style accessibility (p_a) from WCN")
    ap.add_argument("--wt_fasta", required=True)
    ap.add_argument("--wt_pdb", required=True)
    ap.add_argument("--mutants_csv", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--tau", type=float, default=2.0)
    ap.add_argument("--chain", default=None)
    ap.add_argument("--allow_gaps", action="store_true")
    ap.add_argument("--start_idx", type=int, default=None,
                    help="Inclusive WT FASTA start position")
    ap.add_argument("--end_idx", type=int, default=None,
                    help="Inclusive WT FASTA end position")
    args = ap.parse_args()

    wt_seq = read_fasta_one_seq(args.wt_fasta)
    seq_len = len(wt_seq)
    start_idx, end_idx = validate_range(args.start_idx, args.end_idx, seq_len)

    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("wt", args.wt_pdb)

    df_wcn = compute_wcn_allchains(structure, chain_id=args.chain)

    pdb_concat, concat_to_global, _ = concat_chain_sequences(structure, chain_id=args.chain)
    f2p = align_fasta_to_pdbseq(wt_seq, pdb_concat)

    df_wcn["access_raw"] = -df_wcn["wcn"]

    mu = float(df_wcn["access_raw"].mean())
    sd = float(df_wcn["access_raw"].std(ddof=0)) or 1.0
    df_wcn["p_a_pos"] = expit(((df_wcn["access_raw"] - mu) / sd) / float(args.tau))
    df_wcn["p_a_pos"] = np.clip(df_wcn["p_a_pos"], 1e-9, 1 - 1e-9)

    muts = pd.read_csv(args.mutants_csv)
    if "mutant" not in muts.columns:
        raise SystemExit("mutants_csv must contain a 'mutant' column.")

    muts = muts.copy()
    muts["fasta_pos"] = muts["mutant"].apply(parse_pos)

    # Keep only requested WT FASTA range
    muts = muts[
        muts["fasta_pos"].notna() &
        (muts["fasta_pos"] >= start_idx) &
        (muts["fasta_pos"] <= end_idx)
    ].copy()

    def fasta_to_p_a(pos):
        if pos not in f2p:
            return np.nan
        concat_pos = f2p[pos]
        gidx = concat_to_global.get(concat_pos, None)
        if gidx is None:
            return np.nan
        row = df_wcn.loc[df_wcn["global_index"] == gidx]
        if row.empty:
            return np.nan
        return float(row["p_a_pos"].iloc[0])

    muts["p_a"] = muts["fasta_pos"].apply(fasta_to_p_a)

    if args.allow_gaps and not muts.empty:
        muts["p_a"] = muts["p_a"].interpolate(limit_direction="both")
        mval = float(muts["p_a"].mean(skipna=True))
        muts["p_a"] = muts["p_a"].fillna(mval)

    out = muts[["mutant", "p_a"]].dropna().drop_duplicates("mutant")
    out.to_csv(args.out_csv, index=False)
    print(
        f"[accessibility] WT range {start_idx}-{end_idx}; "
        f"wrote {len(out)} rows to {args.out_csv}"
    )


if __name__ == "__main__":
    main()
