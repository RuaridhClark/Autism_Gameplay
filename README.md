# Code Accompanying *Motor Organisation of Social Play in Children with Autism*

This repository contains the data processing and analysis pipeline used in the research paper *"Motor organisation of social play in children with autism"*.

## Repository Overview

The code is organised into modular folders reflecting each stage of the analysis:

### `Track_Swipes/`
- **`identify_swipes.m`** – Processes raw touch-tracking data to detect distinct swipes.  
- **`allocate_swipes.m`** – Organises detected swipes and saves them into individual user files.  

### `Create_adj/`
- **`create_adjs.m`** – Generates adjacency matrices by tracking the zones each swipe passes through.  
  - The adjacency matrix links origin and destination zones during gameplay.  
  - Outputs are stored in the `adjs/` folder.  

### `Metric_creation/`
- **`Sharing_score_analysis.m`** – Computes a "sharing score" for each user by perturbing the adjacency matrix.  
- **`Even_score_analysis.m`** – Evaluates the evenness of eigenvector entries for plate zones, using the standard deviation of entries as the measure.  

### `Results_comparison/`
- **`main_analysis.m`** – Main entry point for running the analysis pipeline.  
  - Loads data, counts swipes, and creates group sets for analysis (e.g., WP, ASD, OND).  
  - Identifies correlations, compares results across groups, and generates plots.  

## Data Notes
- Both pre-trial and trial datasets are supported, but their raw formats differ, requiring different preprocessing steps.  
- Due to data privacy restrictions, **raw swipe data is not included**. Instead, derived network data is provided for reproducibility.  

## Usage
1. Clone this repository.  
2. Add all folders to your MATLAB path.  
3. Run the scripts in the order listed above, or start with `main_analysis.m` in `Results_comparison/` for the full pipeline.  

## License
This project is licensed under the [CC BY 4.0 License](https://creativecommons.org/licenses/by/4.0/).  
