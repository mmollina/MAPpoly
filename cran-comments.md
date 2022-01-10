## Re-submission of MAPpoly package

This is a re-submission of MAPpoly package. This version (0.3.0), contains the following changes

  - Fix minor bugs 
  - Update documentation 
  - Update DESCRIPTION file
  - Added Rccp parallelization to avoid memory overflow in personal computers
  - Added function export_qtlpoly
  - Added function filter_individuals
  - Added several utility functions
  - Added Vignette

Thank you for reviewing our re-submission!

## Test environments
* local R installation (macOS 11.6), R 4.1.2
* local R installation (Ubuntu 20.04), R 4.1.2
* local R installation compiled without long double support (Ubuntu 20.04), R 4.1.1
* Windows Server x64 (on appveyor), R 4.1.2
* Win-builder (4.0.5, 4.1.2, and 4.2.0)

## R CMD check results (on local macOS: R 4.1.2)

0 errors | 0 warnings | 1 note

* installed size is 11.7Mb
  sub-directories of 1Mb or more:
    * R:      2.7Mb
    * data:   3.0Mb
    * doc     4.4Mb

## Downstream dependencies

 There are currently no downstream dependencies for this package
