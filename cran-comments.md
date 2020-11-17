## Re-submission of an archived package

This is a re-submission of an archived package. The package was archived on Nov-16-2020 because we did not correct an error detected by extra checks with long doubles disabled in the required time.

We have corrected and tested this issue in this new version. To test whether the function was working correctly,  we used a Docker container with long doubles disabled with “x86_64 Linux with R-devel” check.

---
Other changes since version 0.2.0

  - Include new test files 
  - Fix minor bugs 
  - Update documentation 
  - Update DESCRIPTION file
  - Add internal functions
  - Remove memuse dependency
  - Added estimation of two-point recombination fraction using genotype probabilities
  - Added 'read_fitpoly' function
  - Updated 'read_vcf' function to import genotype probabilities from PL field 
  - New graphical feature in function 'rf_snp_filt'
  - Added function 'calc_gic': calculates an plot genotypic information content (GIC)
  - Added option to plot large recombination fraction matrix using aggregation of cells of the matrix
  - Added estimation of map distance using the MDSMap package procedure (projects the MDS order onto a single dimension using principal curves)
  

Thank you for reviewing our submission!

## Test environments
* local R installation (macOS 10.15.6), R 4.1.0
* local R installation (Ubuntu 18.04), R 3.6.3
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* Windows Server x64 (on appveyor), R 4.0.2
* Win-builder (3.6.3, 4.0.2, and devel)
* Debian Docker image: "debian-gcc-devel-nold"

## R CMD check results (on local macOS)

0 errors | 0 warnings | 2 note

* New submission; Package was archived on CRAN 

* installed size is 12.2Mb
  sub-directories of 1Mb or more:
    * R:      2.6Mb
    * data:   8.9Mb
    
## Downstream dependencies

 There are currently no downstream dependencies for this package
