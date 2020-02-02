# 2019-nCoV-T-Cell-Vaccine-Candidates

### Prerequisites

Latest version of [**R**](https://www.r-project.org/) (version 3.6.0 or
later) and [**RStudio**](https://rstudio.com/products/rstudio/download/)
(version 1.2.5 or later) installed.

Following **R** packages and their dependencies are required. Upon first
run of `init.R` file, it will check if these are installed and attempt
to install them. Tested to work on **Windows** and **Mac** operating
sytems.

    + here 
    + RColorBrewer 
    + tidyverse 
    + seqinr
    + readxl 
    + readr 
    + glue 
    + reticulate 
    + BALCONY 
    + seqinr 
    + doParallel
    + BiocManager
    + msa 
    + Biostrings
    + doParallel
    + foreach
        
### Organization

-   Main R code is in `main.R` file.
-   Custom functions are in `Functions` folder 
-   All data files are in `Data` folder 
-   Python 2.7 code from [IEDB](tools.iedb.org/population/download/) for computing population coverages is in `Utils` folder

### Main Project File

- Open the `main.R` R script file in RStudio and run. This will run the analysis reproduce the results.
- Outputs will be in `Data` folder.
    
