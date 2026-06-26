# Step 3: TranceptEVE Fitness Prediction

## Overview

This step calculates evolutionary fitness scores for the computational mutation library generated in **Step 2**.

Fitness prediction is performed using the **TranceptEVE** framework, which combines an autoregressive protein language model with evolutionary information from a multiple sequence alignment (MSA). The resulting fitness scores are used as the fitness component (**p<sub>f</sub>**) in the downstream EVEscape-style escape score calculation.

In this pipeline, the **Tranception Large** checkpoint is used together with inference-time retrieval from the repaired MSA generated in Step 2.

---

## Required Inputs

This step requires the following inputs:

| Input | Description |
|--------|-------------|
| Wild-type FASTA | Amino acid sequence of the target antigen |
| Repaired MSA (.a3m) | MSA generated in Step 2 |
| DMS mutation library | Computational DMS library generated in Step 2 |
| Tranception checkpoint | Large pretrained checkpoint downloaded from the official EVEscape repository |

---

## External Dependency

The TranceptEVE model and pretrained checkpoints should be installed following the instructions provided by the official EVEscape repository:

https://github.com/OATML-Markslab/EVEscape

This workflow uses the **Tranception Large** checkpoint:

```
Tranception_Large
```

The model checkpoint is not included in this repository because of its size.

---

## Directory Contents

```
03_trancepteve_fitness/

├── README.md
├── score_tranception_proteingym.py
└── chikv_p62_example_score_tranception_proteingy_0.8.sh
```

---

## Scripts

### score_tranception_proteingym.py

This script performs mutation fitness prediction using the TranceptEVE model.

Main functions include:

- Loading the pretrained Tranception checkpoint
- Reading the repaired MSA generated in Step 2
- Loading the computational DMS mutation library
- Performing inference-time retrieval
- Computing mutation fitness scores
- Writing prediction results for downstream analysis

### Example Shell Script

An example shell script is provided to demonstrate how to run the scoring pipeline.

Run:

```bash
bash chikv_p62_example_score_tranception_proteingy_0.8.sh
```

The example uses:

- Tranception Large checkpoint
- repaired MSA
- computational DMS library
- inference-time retrieval
- retrieval inference weight = **0.8**

---

## Example Command

The shell script executes a command similar to:

```bash
python score_tranception_proteingym.py \
    --checkpoint Tranception_Large \
    --model_framework pytorch \
    --target_seq <wildtype_sequence> \
    --MSA_folder . \
    --MSA_filename <repaired_msa.a3m> \
    --output_scores_folder results \
    --inference_time_retrieval \
    --retrieval_inference_weight 0.8 \
    --batch_size_inference 16 \
    --num_workers 4 \
    --DMS_data_folder . \
    --DMS_file_name <dms_mutants.csv> \
    --MSA_start 1 \
    --MSA_end <protein_length>
```

---

## Retrieval Inference Weight

One of the most important parameters in TranceptEVE is:

```bash
--retrieval_inference_weight
```

When inference-time retrieval is enabled (`--inference_time_retrieval`), TranceptEVE combines information from two complementary sources:

- **Autoregressive transformer predictions**, learned from large protein sequence datasets.
- **MSA-based retrieval**, which incorporates evolutionary information from homologous protein sequences.

The retrieval inference weight controls the relative contribution of the MSA-based retrieval component.

| Retrieval Weight | Interpretation |
|-----------------|----------------|
| 0.2 | Greater reliance on the transformer model with limited contribution from the MSA |
| 0.5 | Balanced contribution from transformer predictions and MSA retrieval |
| 0.8 | Greater reliance on evolutionary information contained in the MSA |

The example shell script included in this repository uses:

```bash
--retrieval_inference_weight 0.8
```

This value is provided **only as an example**.

During development of this pipeline, retrieval weights of **0.2**, **0.5**, and **0.8** were evaluated and compared. The optimal retrieval weight may vary depending on the protein, the quality of the MSA, and the downstream application.

Users are encouraged to experiment with different retrieval weights and select the value that provides the best performance for their own datasets.

---

## Important Parameters

| Parameter | Description |
|-----------|-------------|
| `--checkpoint` | Tranception pretrained checkpoint |
| `--target_seq` | Wild-type protein sequence |
| `--MSA_filename` | Repaired A3M generated in Step 2 |
| `--DMS_file_name` | Computational DMS mutation library |
| `--MSA_start` | Starting residue position of the MSA |
| `--MSA_end` | Ending residue position of the MSA |
| `--batch_size_inference` | Batch size during model inference |
| `--num_workers` | Number of workers for data loading |
| `--inference_time_retrieval` | Enables MSA-guided retrieval during inference |
| `--retrieval_inference_weight` | Weight assigned to the MSA retrieval component |

---

## Output

The prediction results are written to the specified output directory.

Example:

```
results/

└── Tranception_Large_retrieval_0.8_substitutions/
    └── p62_chikv_dms_mutants.csv
```

The output contains mutation-level fitness scores generated by TranceptEVE and serves as the input for **Step 4: Escape Score Calculation**.

---

## Notes

- The repaired A3M generated in Step 2 should be used instead of the original ColabFold MSA.
- The DMS mutation library should contain the `mutated_sequence` column generated in Step 2.
- For full-length proteins, use:

```bash
--MSA_start 1
--MSA_end <protein_length>
```

- GPU acceleration is recommended when scoring large mutation libraries.
- Large pretrained checkpoints should not be committed to Git and should be downloaded separately following the official EVEscape instructions.
