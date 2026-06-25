# Step 4: Escape Score Calculation

## Overview

This step computes the final escape score for each mutation by integrating three complementary components:

* **Evolutionary fitness (p_f)** derived from TranceptEVE
* **Structural accessibility (p_a)** calculated from the wild-type protein structure
* **Physicochemical dissimilarity (p_d)** calculated from amino acid property changes

The final escape score is computed following the EVEscape framework:

[
\text{Escape Score}
===================

\log(p_f)
+
\log(p_a)
+
\log(p_d)
]

An additional score excluding structural accessibility is also generated:

[
\log(p_f)+\log(p_d)
]

This step serves as the core of the antigen escape prediction pipeline.

---

## Required Inputs

The pipeline requires three primary inputs:

| Input              | Description                                          |
| ------------------ | ---------------------------------------------------- |
| Wild-type FASTA    | Amino acid sequence of the target antigen            |
| Wild-type PDB      | Protein structure used for accessibility calculation |
| TranceptEVE output | Mutation fitness predictions generated in Step 3     |

Optionally, users may specify:

* Chain ID
* Sequence region (start and end positions)
* Temperature parameters for each scoring component

---

## Pipeline Components

The workflow is divided into four stages.

### 1. Fitness Component

`fitness_helper.py`

Builds the evolutionary fitness component (**p_f**) from TranceptEVE predictions.

Main tasks include:

* parsing TranceptEVE output
* calculating ΔLLR relative to the wild-type sequence
* converting fitness values into normalized probabilities

Output:

```
fitness_component.csv
```

---

### 2. Accessibility Component

`accessibility_evescape_helper.py`

Calculates residue accessibility (**p_a**) from the wild-type protein structure.

Accessibility is estimated using the EVEscape weighted contact number (WCN) approach.

Main tasks include:

* mapping FASTA positions onto the PDB structure
* computing weighted contact numbers
* converting accessibility values into normalized probabilities

Output:

```
accessibility_component.csv
```

---

### 3. Dissimilarity Component

`dissimilarity_evescape_helper.py`

Calculates mutation dissimilarity (**p_d**) using amino acid physicochemical properties.

Two features are considered:

* Eisenberg–Weiss hydropathy
* Amino acid charge

The combined score is standardized and converted into a probability.

Output:

```
dissimilarity_component.csv
```

---

### 4. Final Score

`evescape_pipeline.py`

The main pipeline script executes the three component calculations and merges their outputs into the final escape score table.

Final output:

```
evescape_scores.csv
```

Output columns include:

| Column         | Description                           |
| -------------- | ------------------------------------- |
| mutant         | Mutation identifier                   |
| p_f            | Fitness component                     |
| p_a            | Accessibility component               |
| p_d            | Dissimilarity component               |
| evescape_score | Composite escape score                |
| p_f+p_d        | Composite score without accessibility |

---

## Running the Pipeline

An example shell script (`example.sh`) is provided to demonstrate how to execute the pipeline with the required inputs and parameters.

Typical usage:

```bash
bash example.sh
```

Alternatively, the pipeline can be executed directly:

```bash
python evescape_pipeline.py \
    --wt_fasta <wildtype.fasta> \
    --wt_pdb <wildtype.pdb> \
    --tr_csv <trancepteve_predictions.csv> \
    --out_dir results/
```

---

## Directory Contents

```
04_escape_scoring/

├── README.md
├── evescape_pipeline.py
├── fitness_helper.py
├── accessibility_evescape_helper.py
├── dissimilarity_evescape_helper.py
└── example.sh
```

---

## Outputs

Running the pipeline produces the following files:

```
results/

├── fitness_component.csv
├── fitness_component.mutants.csv
├── accessibility_component.csv
├── dissimilarity_component.csv
└── evescape_scores.csv
```

The final file (`evescape_scores.csv`) contains the composite escape scores and serves as the input for the downstream analysis and visualization steps.
