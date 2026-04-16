# Morris Water Maze Spatial Learning Simulation Models

*(Place vs Distance Cells)*

---

## Overview

This repository contains simulation and analysis code for comparing:

* Place-cell reinforcement learning model
* Distance-cell model
* Combined place–distance models

Models are evaluated under:

* Fixed platform condition
* Variable platform condition

---

## Repository Structure

```
modelling/             → simulation scripts (generate CSV outputs)
analysis/              → analysis & plotting scripts
plotting_results/      → generated figures
```

---

## Step 1: Run modelling scripts (generate data)

Run scripts in the **`modelling/`** folder:

### Individual models

* `place-cell model.R` → trajectory & cognitive map (place-cell)

* `distance-cell model.R` → trajectory & cognitive map (distance-cell)

* `place-cell model fixed.R` → fixed platform (place-cell)

* `place-cell model variable.R` → variable platform (place-cell)

* `distance-cell model fixed.R` → fixed platform (distance-cell)

* `distance-cell model variable.R` → variable platform (distance-cell)

### Combined model (ratio analysis)

* `place distance cells combined model_mod.R`

  * Runs multiple place–distance ratios
  * Outputs results for both fixed and variable platforms

### Output

These scripts generate:

* Trial-level CSVs (mean ± SE)
* Run-level CSVs (for mixed-effects models)
* Optional `.rds` files (raw PM arrays)

---

## Step 2: Run analysis scripts (generate figures)

Run scripts in the **`analysis/`** folder:

* `comparison_fixed.R`
  → performance plots (fixed platform)

* `comparison_var.R`
  → performance plots (variable platform)

* `combined_comparison_fixed.R`
  → ratio comparison (fixed platform)

* `combined_comparison_var.R`
  → ratio comparison (variable platform)

### Output

These scripts generate:

* Performance plots (PM curves)
* Mixed-effects model tables
* Figures used in the report

---

## Workflow Summary

1. Run modelling scripts → generate CSV outputs
2. Run analysis scripts → generate plots and statistical results

---
