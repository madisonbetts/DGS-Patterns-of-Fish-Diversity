# Reproduce analyses for the Central Interior Highlands Hybridization project \[Gunn et al. 2023]

Follow the steps listed below in the Analyses section to reproduce analyses for this study. Each step below gives a summary of the analysis and directs you to a general code file (e.g., snolh_structure_analysis.Rmd below in Analysis 1) which then works through the analysis step-by-step. This general file will usually point you to other Rmd code, bash shell scripts, or python scripts.

# Project: Central Interior Highlands Hybridization: Spotted, Smallmouth, Neosho, Ouachita, and Little River Bass (SNOLH)

Investigating hybridization and population structure among and within the black basses (genus Micropterus) in the Central Interior Highlands (CIH) ecoregion, including natural and anthropogenic hybridization between Spotted Bass (SPB; Micropterus punctulatus) and species within the Smallmouth Bass species complex (SMBC): Smallmouth Bass (SMB; M. dolomieu), the newly elevated Neosho Bass (NB; M. velox), and two other potentially distinct species in the Ouachita River Basin, the Ouachita Bass (OB; <i>M. cf. dolomieu </i> Ouachita River), and the Little River Bass (LRB; <i>M. cf. dolomieu </i> Little River).</font>

## Abbreviations used in Analysis

Here we give a brief glossary of abbreviations and acronyms used in analyses:

CIH: Central Interior Highlands
SMBC: Smallmouth Bass species complex
SPB: Spotted Bass
SMB: Smallmouth Bass
NB: Neosho Bass
OB: Ouachita Bass
LRB: Little River Bass

## General information on repository structure

This is a publicly visible GitHub repository storing code (and a small amount of data, although we have done our best to avoid uploading large amounts of data due to the limited storage) for Gunn et al. (2023). In the home directory of the repository (SNOLH_Genomics), you will find a README.md file (the source script for this information), the R Project file (SNOLH_Genomics.Rproj), and three different "analysis" directories, each of which corresponds with a specific analysis conducted in our study:

1.  01_filtering_analysis
2.  02_structure_analysis
3.  03_introgress_analysis

Within each analysis directory, you will find an R markdown script (.Rmd) with the name of the analysis, which contains all of the code needed to run the full analysis. Additionally, you will find:

1.  code

The code directory will store all source code, shell scripts, lists of bash commands, and software packages needed for analysis.

Once you have downloaded the repository and located the code directory, you should create two additional sub-directories within each analysis (on the same level as the code directory):

1.  data
2.  figures

The data directory will store all processed data and metadata needed for analysis. The figures folder will contain any raw figures generated in ggplot for each analysis. Ideally, the .Rmd script should have paths set up so that the code reads all data and scripts and generates figures seamlessly.

## Using the code

To reproduce all analyses in Gunn et al. (2023), download this data repository and place in a desired home directory. This may be done on your local machine, but we recommend downloading to a high-performance computing cluster so that all code will run seamlessly in one environment, as long as Rstudio is installed and the GUI can be called on the cluster.

Once all directories are downloaded, create a new sub-directory within the home directory (same level as the seven analysis directories, .Rproj, README.md, etc.) called `/raw_data`. This is where you will store the raw genomic data and associated sample metadata (see <i><b>Data</i></b> section below).

## Data

Raw genotype data and accompanying metadata are available at Zenodo.org: [LINK]

Download these data into to your `/raw_data` directory within the home working directory.

You should have 2 new items in the directory:

1.  snolh_genotype_data.xlsx
2.  snolh_metadata.xlsx

If you have any questions or issues with data and/or code, please don't hesitate to contact me: <jcgunn@uvm.edu>

## Analyses

### Analysis 1: Filtering Analysis

In this analysis, we clean and filter the full genotype data for 487 black bass individuals, which was derived from the diagnostic SNP panel developed by Long et al. (2021) and prepare the data for analysis in Structure and NewHybrids (See Analysis 2). Specifically, we filter the dataset based on three criteria: 1) out poor quality SNP loci (loci that failed to genotype in over 20% of samples); 2) poor quality samples (samples that failed to genotype in over 20% of loci); and 3) potential duplicate samples (samples that are greater than 95% identical across loci).

#### Follow the Code: `01_filtering_analysis/snolh_filtering_analysis.Rmd`

### Analysis 2: Population Structure and Hybrid Assignment

In this analysis, we assess hierarchical population genomic structure among and within Interior Highland species using the diagnostic SNP panel published by Long et al. (2021). We begin with a holistic analysis of population structure among Spotted Bass and all other Interior Highland species (Smallmouth Bass, Neosho Bass, Ouachita Bass, and Little River Bass) and diagnose hybrids between these species using SNPs diagnostic for Spotted Bass. We then exclude detected hybrids and continue with an analysis of population structure and hybridization among Smallmouth Bass, Neosho Bass, Ouachita Bass and Little River Bass. We again exclude any detected hybrids and move on to a final analysis of all Interior Highland species, excluding Spotted Bass and Smallmouth Bass.

#### Follow the Code: `02_structure_analysis/snolh_structure_analysis.Rmd`

### Analysis 3: Introgression Analysis

In this analysis, we further investigate hybridization and introgression within populations that were inferred to contain hybrids based on NEWHYBRIDS analysis (Analysis 2). We use the R package Introgress to regress interspecific heterozygosity on hybrid index for inferred F1, F2, and back-cross individuals at each hierarchical level of hybrid analysis conducted in Analysis 2. With this analysis, we determine whether hybrids are of very recent origin (first or second generation) or if they show a genetic signature of deeper time hybridization. We also infer from this analysis the extent to which non-native alleles have introgressed into the native distribution of each Smallmouth Bass species complex (SMBC) species.

#### Follow the Code: `03_introgress_analysis/snolh_introgress_analysis.Rmd`
