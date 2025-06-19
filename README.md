# Eml4-Alk

Characterizing differences in tumor suppressor function across variants of the EML4-ALK oncogenic fusion

This repository contains code associated with the manuscript

**Diaz-Jimenez A.**, **Shuldiner E. G.**, **Somogyi K.**, **Shih K.**, **Gonzalez O.**, **Najajreh M.**, **Kim S.**, **Akkas F.**, **Murray C. W.**, **Andrejka L.**, **Tsai M. K.**, **Brors B.**, **Hofmann I.**, **Sivakumar S.**, **Sisoudiya S. D.**, **Sokol E. S.**, **Cai H.**, **Petrov D. A.**, **Winslow M. M.**, and **Sotillo R.** EML4-ALK variant-specific genetic interactions shape lung tumorigenesis. Currently available at https://www.biorxiv.org/content/10.1101/2024.08.26.609730v1

# Running the tuba-seq pipeline

For information and code related to running the Tuba-seq pipeline, see https://github.com/eshuldiner/Aging

# Adaptive sampling

To compare equivalent portions of the EML4-ALK V1 and V3 tumor size distribution we use a method which we term “adaptive sampling” in which the same number of tumors per unit of virus delivered are analyzed for each tumor genotype. Specificially we scaled the number of tumors analyzed for each sgRNA i in each cohort j to account for differences in viral titer and the number of mice transduced in each cohort, and then analyzed the largest $N_{i,j}$ tumors per Lenti-sgRNA/Cre vector.

This scaling procedure requires selecting a benchmark sgRNA and a benchmark cohort, and then selecting a defined number of tumors with that benchmark sgRNA from that cohort. The number of tumors sampled for each other sgRNA $i$ in each cohort $j$ ($N_{i,j}$) is then adjusted to take into account the proportions of sgRNAs in the viral pool and differences in the overall viral titer delivered to the young and old cohorts:

```math
N_{i,j} = N_{i=basal, j=basal} * \frac{T_j}{T_{j=basal}} * \frac{p_{i,j}}{p_{i=basal, j=basal}}
```
## Input file

The parameters defining a given analysis are passed through an input file. The parameter file is space-separated, with each line of the file defining an analysis.

By default, the location and name of the project file should be as follows:
<root>/Parameters/<parameter_id>_parameter_file.txt

The following parameters must be provided in this file:

| Parameter Name    | Definition | Datatype |
| -------- | ------- | ------- |
| root  | Root directory  | string |
| project_name | Identifier for the analysis; will prepend output file names   | string |
| analysis    | Type of analysis to run: must be either "basic" (to calculate adaptively sampled percentiles and mean) or "correlation" to    | string |
| focal_dataset | Basal dataset for which a defined number of tumors is specified | string |
| focal_sgID | Basal sgRNA for which a defined number of tumors is specified  | string |
| focal_number | Number of tumors to analyze for the <focal_sgID> in the <focal_dataset>. i.e., $N_{i=basal, j=basal}$ | integer |
| datasets | List of identifiers for datasets to analyze. Indexing must match <datafiles>.  | comma-separated strings |
| datafiles | List of paths to datasets to analyze. Indexing must match < datasets >.  | comma-separated strings |
| info_file | Path to file containing information on samples in <datasets> (including viral titer). See example info file provided (TSG75_samples.txt) | string |
| inerts | Non- or safe-targeting control sgRNAs (used as baseline in calculating the effects of gene inactivation) | comma-separated strings |
| nboot | Number of bootstraps | integer |
| sgids_to_exclude | Optional list of sgIDs to exclude from the analysis | comma-separated strings |
| geneLevel | Optional; if provided statistics will also be calculated at the gene-level (aggregating information across sgRNAs targeting the same gene | boolean |
| rankOn | For correlation analyses; which statistic to use in calculating correlation across contexts. | string |
| incList | For correlation analyses; named list of sgRNAs or genes to include in ranking analysis. Format is name:gene1,gene2...geneN | string |

# Example analysis

An example input file (/Analysis/example.inp) is provided. The files necessary to run this analysis are provided in ./InputFiles. A yaml file to install all necessary dependencies can be found at [https://github.com/eshuldiner/Aging/tree/main/Environment](https://github.com/eshuldiner/Aging/blob/main/Environment/tubaseq.yml).

Usage: 

```
sbatch run_adaptive_array.sh example.inp
```
