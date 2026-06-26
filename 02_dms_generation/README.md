# Step 2: Mutation Library Generation

## Overview

This step prepares the inputs required for TranceptEVE fitness prediction.

Starting from the wild-type antigen sequence and the multiple sequence alignment (MSA) generated in Step 1, the pipeline performs two preprocessing tasks:

1. Generate an **in silico deep mutational scanning (DMS) library** containing all possible single amino acid substitutions.
2. Repair and reformat the ColabFold-generated A3M file so that it is compatible with TranceptEVE.

---

## Required Inputs

| Input                | Description                               |
| -------------------- | ----------------------------------------- |
| Wild-type FASTA      | Amino acid sequence of the target antigen |
| ColabFold MSA (.a3m) | MSA generated in Step 1                   |

---

## Outputs

This step generates three files:

| Output              | Description                                                              |
| ------------------- | ------------------------------------------------------------------------ |
| `*_dms_mutants.csv` | Computational DMS library containing all single amino acid substitutions |
| `*_dms_wt.csv`      | Wild-type sequence in DMS format for TranceptEVE normalization           |
| `*_ready.a3m`       | Cleaned MSA formatted for TranceptEVE                                    |

These files are used directly in **Step 3: TranceptEVE Fitness Prediction**.

---

## Mutation Library Generation

The Python script generates every possible single amino acid substitution across the wild-type sequence.

For a protein of length **L**, the output contains:

* one mutation at every residue position
* substitutions to the remaining 19 standard amino acids

Each mutation record contains:

* mutation label
* residue position
* wild-type residue
* mutant residue
* mutated amino acid sequence

Example:

| mutant | position | wt  | mut |
| ------ | -------: | --- | --- |
| A1C    |        1 | A   | C   |
| A1D    |        1 | A   | D   |
| ...    |      ... | ... | ... |

---

## MSA Preparation

The ColabFold-generated MSA is processed before being used by TranceptEVE.

The preprocessing includes:

* removing lowercase insertion characters
* converting all residues to uppercase
* replacing "." characters with alignment gaps
* removing columns corresponding to gaps in the query sequence
* updating the query header to the format expected by TranceptEVE

The resulting repaired MSA is written as:

```text
*_ready.a3m
```

---

## Running the Pipeline

An example shell script is provided to execute the preprocessing workflow.

Run:

```bash
bash make_dms_and_fix_a3m.sh
```

Alternatively, execute the Python script directly:

```bash
python make_dms_and_fix_a3m.py \
    <wildtype.fasta> \
    <input.a3m> \
    <mutants.csv> \
    <wt.csv> \
    <ready.a3m>
```

---

## Directory Contents

```text
02_dms_generation/

├── README.md
├── make_dms_and_fix_a3m.py
└── make_dms_and_fix_a3m.sh
```

---

## Notes

* The wild-type FASTA sequence should be identical to the sequence used during MSA generation.
* The repaired A3M file should be used for all downstream TranceptEVE calculations.
* The generated DMS library contains only **single amino acid substitutions**, which are the mutations evaluated throughout the remainder of this workflow.
* The wild-type DMS file is generated to facilitate normalization of TranceptEVE fitness scores against the wild-type sequence.
