# Code Accompanying Motor Organisation of Social Play in Children with Autism

This repository contains the data processing and analysis pipeline for the research paper on "Motor organisation of social play in children with autism".

## Directory Structure

### Track_Swipes Folder

- `identify_swipes.m`: Processes raw touch tracking data to identify distinct swipes.
- `allocate_swipes.m`: Organises and saves swipes into a file for each user.

### Create_adj Folder

- `create_adjs.m`: Tracks the zones that swipes pass through to create an adjacency matrix linking the origin and destination zones of swipes during gameplay.

### Metric_creation Folder

- `Sharing_score_analysis.m`: Identifies the sharing score for each user by applying a perturbation to the adjacency matrix.
- `Even_score_analysis.m`: Identifies the evenness of eigenvectors entries for the plate zones by monitoring the standard deviation of said entries.

### Results_comparison Folder

- `main_analysis.m`: Loads data, counts swipes, creates group sets for analysis (e.g. WP, ASD, OND), identifies correlations and compares results, plots results.

## Notes

- The analysis combines pre-trial and trial datasets with variation in raw data formats requiring differences in processing steps.
- Raw swipe data is restricted, with only derived network data included with these files.

## License

This project is licensed under the terms of CC BY 4.0.
