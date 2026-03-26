# Riverscape genetics of non-native Brook Trout to inform native Cutthroat Trout conservation

[https://doi.org/10.5061/dryad.nk98sf84b](https://doi.org/10.5061/dryad.nk98sf84b)

## Description of the data and file structure

This dataset contains the data and sample code necessary to reproduce the riverscape genetic modeling results presented in the paper *"Riverscape genetics of non-native Brook Trout to inform native Cutthroat Trout conservation"*. The study quantifies the effects of riverscape features on trout gene flow by relating pairwise genetic distance between sampling sites to the riverscape features separating them. This analysis applies a statistical framework introduced by White et al. (2020), which represents a stream network as a Spatially Structured Ecological Network (SSEN), a graph consisting of "nodes" (points along a stream network) and "edges" (stream segments connecting nodes). For details on this statistical method, see White et al. (2020).

Empirical genetic and geospatial data are provided along with sample code to implement the riverscape genetic analysis. Genetic data consist of microsatellite genotypes for 757 brook trout (*Salvelinus fontinalis*) collected from 22 sites in the upper Cache la Poudre River watershed, Colorado, USA, with field collections occurring between July and October in 2018 and 2019. Geospatial data include a shapefile representing the study area's ~60 km stream network, comprising 37 stream segments, which serve as "edges" in the SSEN framework. Additionally, measurements of riverscape features for each stream segment are provided as covariates for riverscape genetic modeling, along with a genetic distance matrix representing pairwise genetic differentiation between sampling sites.

### Files and variables

#### File: Stack_et_al_data.zip

**Description:** This is a collection of data files needed to reproduce the riverscape genetic modelling presented in the paper "Riverscape genetics of non-native Brook Trout to inform native Cutthroat Trout conservation". All files must be in the working directory to perform the analysis using the R script “Stack_et_al_analysis_code.r”.

#### Content:

* **stack_et_al_microsatellite_genotypes.rds**: A .csv file containing microsatellite genotypes at 12 loci for 757 individual brook trout collected from 22 sampling sites. 
  * Variables:
    * Sample: Unique identifier for each individual brook trout sample.
    * Site: Unique identifier for each of 22 sampling sites where brook trout were collected. 
    * Remaining columns: paired columns containing diploid genotypes (allele size) for 12 microsatellite loci. Column headers are locus names as described by King et al. 2012, followed by "_1" or "_2" indicating diploid copies. Numeric values for allele sizes represent total fragment length. Missing data are entered as "NA".
* **rwc.fcns.r**: an R source file containing functions necessary to implement the riverscape genetics model using R (see White et al. 2020, Hanks 2018).

- **SSEN_edges.shp, .dbf, .prj, .shx**: geospatial files for viewing a map of the study stream network, separated into 37 stream segments which represent edges in a Spatially Structured Ecological Network (SSEN). This shapefile contains values for riverscape covariates measured on each SSEN edge.
  *     Variables:
    * edgeID: Unique identifier for each edge of the Spatially Structured Ecological Network. 
    * barriers: Binary variable indicating presence (1) or absence (0) of waterfalls on the corresponding edge.
    * damreg: Binary variable indicating whether the corresponding edge is (1) or is not (0) dam-regulated. 
    * gradient: The stream gradient of the corresponding edge (rise/run), standardized by dividing by the maximum value across all edges. 
    * shreve: Shreve's stream order of the corresponding edge, standardized by dividing by the maximum value across all edges.
    * strahler: Strahler's stream order of the corresponding edge, standardized by dividing by the maximum value across all edges.
    * geometry: A MULTILINESTRING geometry column with spatial coordinates in XY format representing a sequence of points defining the stream's path.
- **stack_et_al_covariate_matrices.rds**: an .rds file containing a list of 9 square matrices (38 × 38) representing pairwise covariate values for edges in the SSEN. Matrices are structured such that element *(i, j)* represents the covariate value associated with the edge connecting nodes *i* and *j*. For details on the organization of these covariate matrices to represent an SSEN, see White et al. 2020. 
  * Contents:
    * intercept.cov: Intercept term.
    * gradient.cov: Stream gradient (rise/run) on each edge, standardized by the maximum value across all edges.
    * direction.cov: Binary matrix indicating whether the edge is directed downstream (1) or upstream (0).
    * barriers.cov: Binary matrix indicating the presence (1) or absence (0) of a waterfall on the corresponding edge. 
    * damreg.cov: Binary variable indicating whether the corresponding edge is (1) or is not (0) dam-regulated. 
    * shreve.cov: Shreve's stream order of the corresponding edge, standardized by dividing by the maximum value across all edges.
    * strahler.cov: Strahler's stream order of the corresponding edge, standardized by dividing by the maximum value across all edges.
    * barriers_dir_int.cov: An interaction term between the *barriers* and *direction* covariates.
    * gradient_dir_int.cov: An interaction term between the *gradient *and *direction* covariates.
- **stack_et_al_Fst_matrix.rds:** An .rds file containing a pairwise genetic distance matrix between brook trout sampling sites.  Rows and columns are unique identifiers for brook trout sampling sites, and matrix entries contain estimates of genetic distance between site pairs, measured as Fst (Weir and Cockerham 1984).

## Code/software

#### Script: Stack_et_al_analysis_code.r

This is an R script containing sample code for performing the riverscape genetic analysis presented in the paper "Riverscape genetics of non-native Brook Trout to inform native Cutthroat Trout conservation". All files in the accompanying folder Stack_et_al_data.zip must be in the working directory. The following R packages are required to run the script: "tidyverse" (Wickham et al. 2019), "rwc" (Hanks 2018), and "sf" (Pebesma 2018).

#### References:

Hanks E. M. 2018. rwc: Random Walk Covariance Models. R package version 1.11, [https://CRAN.R-project.org/package=rwc](https://CRAN.R-project.org/package=rwc).

King, T. L., B. A. Lubinski, M. K. Burnham-Curtis, W. Stott, and R. P. Morgan. 2012. Tools for the management and conservation of genetic diversity in Brook Trout (Salvelinus fontinalis): tri- and tetranucleotide microsatellite markers for the assessment of genetic diversity, phylogeography, and historical demographics. Conservation Genetics Resources 4(3):539–543.

Pebesma, E. 2018. Simple Features for R: Standardized Support for Spatial Vector Data. The R Journal 10 (1), 439-446, [https://doi.org/10.32614/RJ-2018-009](https://doi.org/10.32614/RJ-2018-009).

Weir, B. S., and C. C. Cockerham. 1984. Estimating F-Statistics for the Analysis of Population Structure. Evolution 38(6):1358.

White, S. L., E. M. Hanks, and T. Wagner. 2020. A novel quantitative framework for riverscape genetics. Ecological Applications 30(7):e02147.

Wickham H, Averick M, Bryan J, Chang W, McGowan LD, François R, Grolemund G, Hayes A,Henry L, Hester J, Kuhn M, Pedersen TL, Miller E, Bache SM, Müller K, Ooms J, Robinson D, Seidel DP, Spinu V, Takahashi K, Vaughan D, Wilke C, Woo K, Yutani H (2019).“Welcome to the tidyverse.” \_Journal of Open Source Software, 4(43), 1686. doi:10.21105/joss.01686. [https://doi.org/10.21105/joss.01686](https://doi.org/10.21105/joss.01686).