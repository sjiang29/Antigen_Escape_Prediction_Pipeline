# Antigen Escape Prediction Pipeline

This repository provides an end-to-end workflow for predicting potential antibody escape mutations from antigen sequences.

The pipeline integrates:

* Evolutionary information from multiple sequence alignments (MSAs)
* Fitness prediction using TranceptEVE
* Structure-based residue accessibility
* Amino acid physicochemical dissimilarity
* Downstream statistical analysis
* Structural visualization in PyMOL

## Workflow

```text
Antigen Sequence
      │
      ▼
01. MSA Generation (ColabFold/MMseqs2)
      │
      ▼
02. Mutation Library Generation
      │
      ▼
03. TranceptEVE Fitness Prediction
      │
      ▼
04. Escape Scoring
      │
      ▼
05. Data Analysis
      │
      ▼
06. Structural Visualization
```

## Repository Organization

| Step | Directory                | Description                           |
| ---- | ------------------------ | ------------------------------------- |
| 01   | `01_msa_generation`      | Generate MSAs using ColabFold         |
| 02   | `02_dms_generation`      | Generate mutation libraries           |
| 03   | `03_trancepteve_fitness` | Calculate mutational fitness scores   |
| 04   | `04_escape_scoring`      | Calculate composite escape scores     |
| 05   | `05_analysis`            | Statistical analysis and benchmarking |
| 06   | `06_visualization`       | PyMOL visualization scripts           |

Each directory contains:

* Detailed README
* Example inputs
* Example outputs
* Source code
* Usage instructions

## Getting Started

Begin with:

```text
01_msa_generation/README.md
```

and follow the workflow sequentially through each stage.

## License

[Choose appropriate license]

## Citation

If you use this repository, please cite the associated publication (when available).
