#!/usr/bin/env python3
"""
evescape_pipeline.py
--------------------
Orchestrator that builds the three EVEscape components (p_f, p_a, p_d) from:
  - WT FASTA      (wt.fasta)
  - WT PDB        (wt.pdb)
  - TranceptEVE   (trancepteve.csv with mutated_sequence & avg_score columns)

Outputs:
  mutant, p_f, p_a, p_d, evescape_score, p_f+p_d

Optional region restriction:
  --start_idx / --end_idx are inclusive WT FASTA positions.
  They are forwarded to accessibility and dissimilarity helpers so the final
  merged output only contains mutants in the requested WT sequence region.
"""
import argparse
import subprocess
from pathlib import Path
import sys
import pandas as pd
import numpy as np


def run(cmd):
    print("[run]", " ".join(cmd))
    ret = subprocess.run(cmd)
    if ret.returncode != 0:
        sys.exit(f"Command failed with code {ret.returncode}: {' '.join(cmd)}")


def main():
    ap = argparse.ArgumentParser(description="Build EVEscape scores from three inputs")
    ap.add_argument("--wt_fasta", required=True)
    ap.add_argument("--wt_pdb", required=True)
    ap.add_argument("--tr_csv", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--sep", default=None, help="Separator for TranceptEVE CSV (auto if None)")
    ap.add_argument("--tau_f", type=float, default=2.0)
    ap.add_argument("--tau_a", type=float, default=2.0)
    ap.add_argument("--tau_d", type=float, default=2.0)
    ap.add_argument("--chain", default=None, help="Chain ID for PDB (optional)")
    ap.add_argument("--start_idx", type=int, default=None,
                    help="Inclusive WT FASTA start position")
    ap.add_argument("--end_idx", type=int, default=None,
                    help="Inclusive WT FASTA end position")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # 1) FITNESS
    fitness_csv = str(out_dir / "fitness_component.csv")
    mutants_csv = str(out_dir / "fitness_component.mutants.csv")
    cmd_fitness = [
        sys.executable, str(Path(__file__).with_name("fitness_helper.py")),
        "--wt_fasta", args.wt_fasta,
        "--tr_csv", args.tr_csv,
        "--out_csv", fitness_csv,
        "--tau", str(args.tau_f),
    ]
    if args.sep:
        cmd_fitness += ["--sep", args.sep]
    if args.start_idx is not None:
        cmd_fitness += ["--start_idx", str(args.start_idx)]
    if args.end_idx is not None:
        cmd_fitness += ["--end_idx", str(args.end_idx)]
    run(cmd_fitness)

    # 2) ACCESSIBILITY
    access_csv = str(out_dir / "accessibility_component.csv")
    cmd_access = [
        sys.executable, str(Path(__file__).with_name("accessibility_evescape_helper.py")),
        "--wt_pdb", args.wt_pdb,
        "--wt_fasta", args.wt_fasta,
        "--mutants_csv", mutants_csv,
        "--out_csv", access_csv,
        "--tau", str(args.tau_a),
    ]
    if args.chain:
        cmd_access += ["--chain", args.chain]
    if args.start_idx is not None:
        cmd_access += ["--start_idx", str(args.start_idx)]
    if args.end_idx is not None:
        cmd_access += ["--end_idx", str(args.end_idx)]
    run(cmd_access)

    # 3) DISSIMILARITY
    dissim_csv = str(out_dir / "dissimilarity_component.csv")
    cmd_dissim = [
        sys.executable, str(Path(__file__).with_name("dissimilarity_evescape_helper.py")),
        "--mutants_csv", mutants_csv,
        "--out_csv", dissim_csv,
        "--tau", str(args.tau_d),
    ]
    if args.start_idx is not None:
        cmd_dissim += ["--start_idx", str(args.start_idx)]
    if args.end_idx is not None:
        cmd_dissim += ["--end_idx", str(args.end_idx)]
    run(cmd_dissim)

    # 4) COMBINE
    f = pd.read_csv(fitness_csv)
    a = pd.read_csv(access_csv)
    d = pd.read_csv(dissim_csv)

    if "mutant" not in f.columns or "p_f" not in f.columns:
        sys.exit("fitness_component.csv must contain columns: mutant, p_f")
    if "mutant" not in a.columns or "p_a" not in a.columns:
        sys.exit("accessibility_component.csv must contain columns: mutant, p_a")
    if "mutant" not in d.columns or "p_d" not in d.columns:
        sys.exit("dissimilarity_component.csv must contain columns: mutant, p_d")

    out = (
        f[["mutant", "p_f"]]
        .merge(a[["mutant", "p_a"]], on="mutant", how="inner")
        .merge(d[["mutant", "p_d"]], on="mutant", how="inner")
    )

    out["evescape_score"] = np.log(out["p_f"]) + np.log(out["p_a"]) + np.log(out["p_d"])
    out["p_f+p_d"] = np.log(out["p_f"]) + np.log(out["p_d"])

    final_csv = str(out_dir / "evescape_scores.csv")
    out.to_csv(final_csv, index=False)
    print(f"[final][done] {len(out)} mutants written to {final_csv}")
    print("Columns: mutant, p_f, p_a, p_d, evescape_score, p_f+p_d")


if __name__ == "__main__":
    main()
