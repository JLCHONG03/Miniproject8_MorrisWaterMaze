Morris Water Maze Spatial Learning Simulation Models (Place vs Distance Cells) 

This repository contains simulation and analysis code for comparing:
•	Place-cell reinforcement learning model
•	Distance-cell model
•	Combined place–distance models

Models are evaluated under:
•	Fixed platform condition
•	Variable platform condition
________________________________________

Repository Structure
modelling/   → simulation scripts (generate CSV outputs)
analysis/    → analysis & plotting scripts (use CSV outputs)
plotting_results/   → generated visualizations plot

________________________________________

Step 1: Run modelling scripts (generate data)
Run the following scripts in the modelling/ folder to generate CSV files:
Individual models
•	place-cell model.R → trajectory and cognitive map (place-cell)
•	distance-cell model.R → trajectory and cognitive map (distance-cell)
•	place-cell model fixed.R → fixed platform (place-cell)
•	place-cell model variable.R → variable platform (place-cell)
•	distance-cell model fixed.R → fixed platform (distance-cell)
•	distance-cell model variable.R → variable platform (distance-cell)
Combined models (ratio analysis)
•	place distance cells combined model_mod.R
→ generates results for multiple place–distance ratios
→ outputs CSV files for both fixed and variable platform conditions
Output
These scripts will generate:
•	trial-level summary CSVs (mean ± SE)
•	run-level CSVs (for mixed-effects models)
•	optionally .rds files (raw PM arrays)

________________________________________

Step 2: Run analysis scripts (generate figures)
After generating CSV files, run scripts in the analysis/ folder:
•	comparison_fixed.R
→ plots performance measures for fixed platform (place vs distance)
•	comparison_var.R
→ plots performance measures for variable platform (place vs distance)
•	combined_comparison_fixed.R
→ plots performance across different ratios (fixed platform)
•	combined_comparison_var.R
→ plots performance across different ratios (variable platform)
Output
These scripts generate:
•	Performance plots (PM curves)
•	Mixed-effects model result tables
•	Figures used in the report
________________________________________

Workflow Summary
1.	Run modelling scripts → generate CSV outputs
2.	Run analysis scripts → generate plots and statistical results
________________________________________
