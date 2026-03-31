# Data from: The role of ecology in allopatric speciation of darters in the Central Highlands, USA

Dataset DOI: [10.5061/dryad.stqjq2cfv](https://doi.org/10.5061/dryad.stqjq2cfv)

## Description of the data and file structure

The dataset includes data and scripts to reproduce the analyses in Stokes et al., 2025 Evolution. In this project, we analyze double-digest restriction-site associated DNA sequence data (archived at the NCBI archive, BioProject Accession Number PRJNA1177478). We perform phylogenetic and population genomic analyses. We then analyze abiotic variables to test for evidence of ecological niche divergence. We also analyze morphological and meristic trait data and gut contents from museum specimens. There are diverse data formats and software used throughout these analyses.

### Files and variables

#### File: Tree_Files.zip

**Description:** This zip folder contains two files. The trees can be viewed in FigTree.

1. "IQ-Tree.tre" - the maximum-likelihood tree from the concatenated ddRAD dataset.
2. "Snapper.tree" - the time-calibrated species tree inferred using SVDQuartets and the snapper module in Beast.

#### File: Meristics_and_Morphology.zip

**Description:** Datasets for analyzing meristic trait data and snout morphology.

1. "Meristic_Analysis.R" --> script for analyzing meristic trait data including PCA, Mahalanobis Distance calculation and Linear Discriminant Analysis.
   1. A file of meristic trait data for each putative species (9 files): "Meristics_P_speciesname.xslx". This tabular datasets includes the following variables:

      | Variable   | Description                                                      | Unit                                                                                                                                                        |
      | :--------- | :--------------------------------------------------------------- | :---------------------------------------------------------------------------------------------------------------------------------------------------------- |
      | Catalog    | Museum Code and Accession Number                                 | Museum\_Number                                                                                                                                              |
      | Drainage   | The river drainage where the individual fish was collected from. | River basins are all in the southeastern United States.                                                                                                     |
      | Individual | Denotes the exact individual used for the measurement.           | Number denotes the exact individual. Associated tag number is noted if available. Also noted is if the individual is a paratype, holotype, or paratopotype. |
      | Sex        | Sex                                                              | F = female; M = male                                                                                                                                        |
      | SL         | Standard Length                                                  | centimeters                                                                                                                                                 |
      | LL         | Lateral Line                                                     | Number of Scales                                                                                                                                            |
      | PoreLL     | Pored Lateral Line                                               | Number of Scales                                                                                                                                            |
      | AbLL       | Scale rows above the lateral line                                | Number of rows                                                                                                                                              |
      | BlwLL      | Scale rows below the lateral line                                | Number of rows                                                                                                                                              |
      | Trans      | Transverse scale rows                                            | Number of rows                                                                                                                                              |
      | CD         | Scales around the caudal peduncle                                | Number of scales                                                                                                                                            |
      | D1         | Dorsal Fin Spines                                                | Number of Spines                                                                                                                                            |
      | D2         | Dorsal Fin Rays                                                  | Number of Rays                                                                                                                                              |
      | P1         | Pectoral Fin Rays                                                | Number of Rays                                                                                                                                              |
      | A1         | Anal Fin spines                                                  | Number of Spines                                                                                                                                            |
      | A2         | Anal Fin rays                                                    | Number of rays                                                                                                                                              |

      <br />
2. "Snout_Data.xlsx", snout length, standard length, and head length measurements. This tabular datasets includes the following variables (if not previously noted above).

| Variable | Description                                                   | Unit                                                                                               |
| :------- | :------------------------------------------------------------ | :------------------------------------------------------------------------------------------------- |
| Species  | Lineage identity as described in the accompanying manuscript. | brucethompsoni = *Percina brucethompsoni;* nasuta\_Arkansas = *Percina cf. nasuta* Arkansas (etc.) |
| SnL/SL   | Ratio of snout length to standard length                      | Ratio (cm/cm)                                                                                      |
| Snl/HL   | Ratio of snout length to head length                          | Ratio (cm/cm)                                                                                      |

<br />

#### File: Diet_Analysis.zip

**Description:** Two files including the raw data for diet and an R script for processing.

1. "Diet_Analysis.R": file for NMDS and PERMANOVA analysis in R.
2. "Percina_swainiadiet2.csv" --> raw gut content data. This tabular datasets includes the following variables (if not previously noted above).

   | Variable                                                    | Description                                                                                                        | Unit                                                                                             |
   | :---------------------------------------------------------- | :----------------------------------------------------------------------------------------------------------------- | :----------------------------------------------------------------------------------------------- |
   | Code                                                        | code used in R script for species identity of each individual                                                      | Corresponds to lineages of *Percina phoxocephala, Percina nasuta,* and *Percina brucethompsoni.* |
   | Month                                                       | Month collected                                                                                                    | 1-12 (1 = January, etc).                                                                         |
   | Stomach Label Number                                        | Specific individual examined                                                                                       | n/a                                                                                              |
   | N                                                           | Denotes whether the individual had a full or empty stomach                                                         | 1 = full; 0 = empty                                                                              |
   | Standard Length                                             | n/a                                                                                                                | cm                                                                                               |
   | ODANATA, Aeshnidae, Coenagrionidae, and remaining variables | The number of each prey item observed, the column names represent the prey item and the numbers represent a count. | Count                                                                                            |

#### File: Ecological_Niche_Models.zip

**Description:** Nine files and scripts to reproduce the ecological niche modeling. Our locality data was combined with locality information from Fishnet2, we snapped these localities to the National Hydrography Dataset (NHD) in ArcGIS Pro to assign each of the localities a unique flowline ID (COMID). The first six files are used to query the EPA's streamcat database and compile/clean the results. Files 7-9 can be used for analysis. If you would simply like to reproduce the analysis using the pre-compiled data start with files 7-9.

1. "ddRADstudyspeciesNHD.csv" - localiites from our study snapped to NHD flowlines. This tabular datasets includes the following variables (if not previously noted above).

| Variable             | Description                                                       | Unit              |
| :------------------- | :---------------------------------------------------------------- | :---------------- |
| Clone                | Clone code for tissues used in analyses.                          | n/a               |
| Name Code            | Lineage assigned to each individual                               | n/a               |
| Locality Description | Description of where the collection was made                      | text              |
| Latitude             | n/a                                                               | Decimal degrees   |
| Longitude            | n/a                                                               | Decimal degrees   |
| COMID                | Unique identifier for the NHD flowline the locality is closest to | n/a               |
| LENGTHKM             | length of NHD flowline                                            | km                |
| AreaSqKM             | area of catchment                                                 | km^2^             |
| TotDASqKM            | drainage area at flowline                                         | km^2^             |
| MAXELEVSMO           | maximum flowline elevation                                        | cm                |
| MINELEVSMO           | minimum flowline elevation                                        | cm                |
| SLOPE                | slope of flowline                                                 | meters/meters     |
| QE\_MA               | Mean Annual Flow from gage adjustment                             | ft^3^/s           |
| HUC8                 | the HUC8 drainage basin code associated with each locality        | unique identifier |

1. "SwainiaNHDFlowline.csv" - localities from Fishnet2 snapped to NHD flowlines. Variables are as above.
2. "randombackground30percent.csv" - randomly sampled COMID's from the background range of each species (sampled using spatially uniform distribution). Variables as are above.
3. "rockstats3.csv" - results of catchment level rock-type percentages. Compiled from the State Geologic Map Compilation (SGMC) and spatial analyses done in ArcGIS Pro. Variables include COMID associated with each catchment (as above), the rock-type, and percentage area in catchment underlain by that rock-type.
4. "StreamCataVariableList2.csv" - list of streamcat variables to query. The names are used to query the database, the description describes the variable, and the UseMe variable indicates whether it is queried in our analyses (1 = yes, 0 = no).
5. "getStreamCataData.R" - use this script to query the streamcat database and prepare/clean data
6. "ENM_Data.csv" output from the R script in #6. Variables are as follows:

   | Variable                                            | Description                                                | Units         |
   | :-------------------------------------------------- | :--------------------------------------------------------- | :------------ |
   | Clastic, Shale, Igneous, Chert, Detrital, Carbonate | Percentage of catchment underlain by the rock-type listed  | %             |
   | PRECIP8110CAT                                       | Mean annual precipitation in catchment (1981-2011)         | mm            |
   | TMAX8110CAT                                         | Mean maximum annual precipitation in catchment (1981-2011) | degrees C     |
   | SANDCAT                                             | Mean % sand in soils in catchment                          | %             |
   | Forest                                              | Percentage of catchment forested                           | %             |
   | Slope                                               | Log10(slope)                                               | meters/meters |
   | RCKDEPCAT                                           | Mean depth to bedrock in catchment                         | cm            |
   | QE\_MA                                              | log10(mean annual discharge)                               | ft^3^/s       |
7. "ENM_MAXENT.R" - R script to reproduce species distribution modelling
8. "ENM_PCA.R" - R script to conduct PCA

#### File: Genomic_Analyses.zip

**Description:** Files and scripts to reproduce genomic analyses.

1. Eight files for G-PHoCS analysis (data and scripts).
2. "Swainia_snapper.xml" - used in snapper module of BEAST for estimating divergence times
3. Two files for reproducing SVDQuartet analysis
4. "Swainia_TreeMix.hdf5" data to reproduce TreeMix analysis.
5. Three control files for BPP analysis (gdi estimation).
6. An example ipyrad parameter file.

\**Subfolder IBD: **Files and scripts for isolation-by-distance analysis.

1. "IBD_data.vcf" - filtered loci used for genetic distance estimation. Data was filtered in VCFTools to retain only biallelic sites that contained no more than 80% missing data. We then thinned our data by sampling once per 10,000 base pairs resulting in a dataset with 8,848 loci.
2. "streamwise_distance.m" MATLAB script for finding pairwise streamwise distance between localities. This accepts a digital elevation model input.
3. "streamwise_distance Fxn.m" MATLAB function for calculating pairwise streamwise distance between localities using flow-routing tools from TopoToolbox.
4. "vcf_pops.csv" location and species identity of each locality. Variables are as above.
5. "streamwisedist2.csv" results of MATLAB streamwise distance calculation.
6. "IBD.R" script for testing for evidence of IBD within and between putative species.

## Code/software

**Genomic Analyses**

ipyrad ([https://ipyrad.readthedocs.io/en/master/](https://ipyrad.readthedocs.io/en/master/))

IQ-TREE ([https://iqtree.github.io/](https://iqtree.github.io/))

R version 4.4.1

G-PhoCS

VCFTools

Beast2

SVDQuartets

BPP

**Isolation-by-distance**

MATLAB

TopoToolbox2 ([https://topotoolbox.wordpress.com/download/](https://topotoolbox.wordpress.com/download/))

## Access information

Other publicly accessible locations of the data:

* DNA sequence data is available at the NCBI archive - BioProject Accession Number PRJNA1177478.

Data was derived from the following sources:

* EPA StreamCat ([https://www.epa.gov/national-aquatic-resource-surveys/streamcat-dataset](https://www.epa.gov/national-aquatic-resource-surveys/streamcat-dataset))
* Fishnet2 ([https://www.fishnet2.net/](https://www.fishnet2.net/))
* State Geologic Map Compilation ([https://www.sciencebase.gov/catalog/item/5888bf4fe4b05ccb964bab9d](https://www.sciencebase.gov/catalog/item/5888bf4fe4b05ccb964bab9d))
* National Hydrography Dataset v2 ([https://www.usgs.gov/national-hydrography/national-hydrography-dataset](https://www.usgs.gov/national-hydrography/national-hydrography-dataset))
* Hydrosheds hydrologically conditioned North America digital elevation model (3-arc second) ([https://www.hydrosheds.org/products](https://www.hydrosheds.org/products))

