## Re-submission of archived package

This is a re-submission of an archived package. The package was archived on Nov-16-2020 because we did not correct an error detected by extra checks with long doubles disabled in the required time.

We have corrected and tested this issue in this new version. To test whether the function was working correctly, we compiled an R version 4.0.3 without long double support in an Ubuntu 18.04 machine. We confirmed the lack of long double support by calling the function capabilities() in R.  

---
Other changes since version 0.2.0

  - Include new test files 
  - Fix minor bugs 
  - Update documentation 
  - Update DESCRIPTION file
  - Add internal functions
  - Remove 'memuse' dependency
  - Added new features in functions 'est_pairwise_rf', 'reest_rf', 'rf_snp_filt' and 
    'plot.mappoly.rf.matrix' (see NEWS.md for details)
  - Added 'read_fitpoly' function

Thank you for reviewing our submission!

## Test environments
* local R installation (macOS 10.15.6), R 4.1.0
* local R installation (Ubuntu 18.04), R 3.6.3
* local R installation compiled without long double support (Ubuntu 18.04), R 4.0.2
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* Windows Server x64 (on appveyor), R 4.0.2
* Win-builder (3.6.3, 4.0.2, and devel)

## R CMD check results (on local macOS)

0 errors | 0 warnings | 2 note

* New submission; Package was archived on CRAN 

* installed size is 12.2Mb
  sub-directories of 1Mb or more:
    * R:      2.6Mb
    * data:   8.9Mb
    
## Downstream dependencies

 There are currently no downstream dependencies for this package
