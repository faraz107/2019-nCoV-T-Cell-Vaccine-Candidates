# 2019-nCoV-T-Cell-Vaccine-Candidates

### Prerequisites

Latest version of [**R**](https://www.r-project.org/) (version 3.6.0 or
later) and [**RStudio**](https://rstudio.com/products/rstudio/download/)
(version 1.2.5 or later) installed.

Set the working directory to the downloaded repository folder in local system.

Following **R** packages and their dependencies are required. 

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
-   Scripts for preparation of raw genomic sequences is in `preparation.R`
-   Custom functions are in `Functions` folder 
-   All data files are in `Data` folder 
-   Python 2.7 code from [IEDB](http://tools.iedb.org/population/download/) for computing population coverages is in `Utils` folder
-   MAFFT v7 source code for Linux operating system downloaded from [here]( https://mafft.cbrc.jp/alignment/software/linux.html)

### Main Project File

- Open the `main.R` R script file in RStudio and run. This will run the analysis reproduce the results.
- Outputs will be in `Data` folder.
- `Preparation.R` will work only with Linux MAFFT software installed and available on the system.
- *Population coverage code* works only with an available Conda Python 2.7 environment (named as "py27") on the system.
- If the said Conda environment is not avaialable then the population coverages can be computed online using [IEDB Analysis Resource](http://tools.iedb.org/population/). 
    - Select the "Population" (*China* or *World*).
    - Select "Calculation option" (*Class I and II combined*) 
    - Proide one of the following files (located within `Data` folder) to the option "Enter epitope / MHC restriction data in the form below or select a file":
        - `Tcell_epitopes_China` 
            (for computing coverages for set of MHC alleles identified to maximize population coverage in China.)
        - `Tcell_epitopes_World` 
            (for computing coverages for set of MHC alleles identified to maximize global population coverage.)
 
 #### [Detailed steps for downloading raw data](https://faraz107.gitlab.io/wuhanvirusanalysis/)
 
 #### [Processed MSAs](https://faraz107.gitlab.io/wuhanvirusanalysis/)
 
 #### [Detailed steps for population coverage computation through IEDB online resource](https://faraz107.gitlab.io/wuhanvirusanalysis/)
 
 ### Troubleshooting
 
 For any questions or comments related to the code and data, please send email to: faraz107@gmail.com.
