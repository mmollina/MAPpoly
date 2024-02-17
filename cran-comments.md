## Re-submission of MAPpoly package

This is a re-submission of MAPpoly package. This version (0.4.0), contains the following changes

  - Update DESCRIPTION file
  - Functions to build maps in individual parents
  - Functions to merge individual maps
  - Imputation functions based on map
  - Functions to edit order interactively
  - Fix minor bugs 

Thank you for reviewing our re-submission!

## Test environments
* local R installation (macOS Somona 14.3.1), R 4.3.2
* local R installation (Ubuntu 22.04), R 4.3.2
* local R installation (Ubuntu 22.04), R devel with  
  - Intel(R) oneAPI DPC++/C++ Compiler 2024.0.2 (2024.0.2.20231213)
  - ifx (IFX) 2024.0.2 20231213
* Win-builder (release and devel)

## R CMD check results 

0 errors | 0 warnings | 2 notes

 - on local macOS: R 4.3.2
   * installed size is 8.2Mb
     sub-directories of 1Mb or more:
       * R:      2.8Mb
       * data:   3.0Mb
   * GNU make is a SystemRequirements.      
       
## Downstream dependencies

We checked 3 reverse dependencies (3 from CRAN + 0 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
