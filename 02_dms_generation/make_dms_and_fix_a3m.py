#!/usr/bin/env python3
import sys
import csv
import os

AA_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

def read_fasta_sequence(fasta_path):
    """Return the first sequence (concatenated, uppercase) from a FASTA."""
    seq = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            s = line.strip()
            if s:
                seq.append(s)
    s = "".join(seq).upper()
    if not s:
        raise ValueError(f"No sequence found in FASTA: {fasta_path}")
    bad = [c for c in s if c not in AA_ALPHABET]
    if bad:
        raise ValueError(f"WT contains non-standard letters {sorted(set(bad))}. Allowed: {AA_ALPHABET}")
    return s

def load_a3m_as_records(a3m_path):
    """Load A3M into list of (header, seq_string) records, preserving order."""
    records = []
    header = None
    seq_chunks = []
    with open(a3m_path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(seq_chunks)))
                header = line
                seq_chunks = []
            else:
                seq_chunks.append(line)
    if header is not None:
        records.append((header, "".join(seq_chunks)))
    if not records or not records[0][0].startswith(">"):
        raise ValueError("A3M must start with a '>' header line.")
    return records

def clean_a3m_letters(seq):
    """Uppercase, drop lowercase (insertions), replace '.' with '-'."""
    out = []
    for ch in seq:
        if "a" <= ch <= "z":
            continue
        if ch == ".":
            out.append("-")
        else:
            out.append(ch.upper())
    return "".join(out)

def mask_columns_by_query(records):
    """
    Remove alignment columns where the (cleaned) query has '-' in that column.
    Return new records (same headers) and the masked query string.
    """
    headers = [h for h,_ in records]
    seqs = [clean_a3m_letters(s) for _, s in records]

    query = seqs[0]
    keep_idx = [i for i, c in enumerate(query) if c != "-"]
    if not keep_idx:
        raise ValueError("Query becomes empty after masking; check your A3M formatting.")

    def apply_mask(s):
        return "".join(s[i] for i in keep_idx)

    masked = [apply_mask(s) for s in seqs]
    new_records = list(zip(headers, masked))
    query_masked = masked[0]
    return new_records, query_masked

def write_a3m(records, out_path, force_first_header=None):
    """Write records to A3M; sequence lines unwrapped (one line per record)."""
    with open(out_path, "w") as f:
        for idx, (h, s) in enumerate(records):
            if idx == 0 and force_first_header is not None:
                h = force_first_header
            f.write(h + "\n")
            f.write(s + "\n")

def make_all_single_mutants(wt):
    """Yield dict rows for every single-AA substitution across the WT."""
    L = len(wt)
    for i, wt_aa in enumerate(wt, start=1):  # 1-based
        for mut_aa in AA_ALPHABET:
            if mut_aa == wt_aa:
                continue
            mut_seq = wt[:i-1] + mut_aa + wt[i:]
            yield {
                "mutant": f"{wt_aa}{i}{mut_aa}",
                "position": i,
                "wt": wt_aa,
                "mut": mut_aa,
                "mutated_sequence": mut_seq
            }

def dedup_by_mutated_sequence(rows):
    """Remove duplicates by 'mutated_sequence' while keeping first occurrence."""
    seen = set()
    out = []
    for r in rows:
        key = r["mutated_sequence"]
        if key in seen:
            continue
        seen.add(key)
        out.append(r)
    return out

def write_mutants_csv(rows, out_path):
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["mutant","position","wt","mut","mutated_sequence"])
        writer.writeheader()
        for r in rows:
            writer.writerow(r)

def write_wt_csv(wt, out_path):
    """
    Write a 1-row DMS CSV for WT that includes a dummy 'mutant' column.
    Use midpoint i with no-op AA->same AA, so Tranception's windowing works.
    """
    L = len(wt)
    mid = L // 2 + 1  # 1-based midpoint
    aa = wt[mid - 1]
    row = {
        "mutant": f"{aa}{mid}{aa}",
        "position": mid,
        "wt": aa,
        "mut": aa,
        "mutated_sequence": wt
    }
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["mutant","position","wt","mut","mutated_sequence"])
        writer.writeheader()
        writer.writerow(row)

def main():
    if len(sys.argv) != 6:
        print("Usage: python make_dms_and_fix_a3m.py <wt.fasta> <msa.a3m> <out_mutants.csv> <out_wt.csv> <out_repaired_a3m>")
        sys.exit(1)

    wt_fasta, msa_a3m, out_mut_csv, out_wt_csv, out_repaired_a3m = sys.argv[1:6]

    # 1) Read WT
    wt = read_fasta_sequence(wt_fasta)
    L = len(wt)

    # 2) Load & clean MSA, mask columns where query has '-' so query becomes gap-free
    records = load_a3m_as_records(msa_a3m)
    cleaned_records = [(h, clean_a3m_letters(s)) for (h, s) in records]
    masked_records, masked_query = mask_columns_by_query(cleaned_records)

    # 3) Sanity checks vs WT
    if len(masked_query) != L:
        print(f"[WARN] Masked query length ({len(masked_query)}) != WT length ({L}). "
              f"We'll still generate mutants from the WT FASTA. "
              "Ensure you pass --MSA_start 1 --MSA_end WT_len when scoring.")
    if masked_query != wt:
        print("[WARN] Masked query != WT sequence. We'll keep WT from FASTA as truth for DMS generation.\n"
              "       Consider aligning your A3M query to exactly match WT (uppercase, no gaps).")

    # 4) Force first header to EVE/Tranception-style span
    forced_header = f">FOCUS/1-{L}"
    masked_records[0] = (forced_header, masked_records[0][1])

    # 5) Write repaired A3M
    os.makedirs(os.path.dirname(out_repaired_a3m) or ".", exist_ok=True)
    write_a3m(masked_records, out_repaired_a3m, force_first_header=forced_header)
    print(f"✅ Repaired A3M written to: {out_repaired_a3m}")

    # 6) Build computational DMS of all single substitutions from WT
    rows = list(make_all_single_mutants(wt))
    rows = dedup_by_mutated_sequence(rows)
    write_mutants_csv(rows, out_mut_csv)
    print(f"✅ Mutants DMS written to: {out_mut_csv}  (rows: {len(rows)})")

    # 7) Write WT-only CSV WITH a dummy 'mutant' column
    write_wt_csv(wt, out_wt_csv)
    print(f"✅ WT DMS written to: {out_wt_csv}")

    # 8) Final tips
    print("\nTips:")
    print("  • Use the repaired A3M for retrieval:  --MSA_filename", out_repaired_a3m)
    print(f"  • Set window to full WT:               --MSA_start 1 --MSA_end {L}")
    print("  • Run mutants and WT separately, then normalize ΔLL vs WT.")

if __name__ == "__main__":
    main()
