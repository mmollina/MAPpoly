## Re-submission of MAPpoly package

This is a re-submission of MAPpoly package. In this version (0.2.3), suggested packages are used conditionally, following §1.1.3.1 of 'Writing R Extensions'. We also added the flag LazyDataCompression: xz in DESCRIPTION file, following §1.1.6 of 'Writing R Extensions'. 

---
Other changes since version 0.2.1

  - Fix minor bugs 
  - Update documentation 
  - Update DESCRIPTION file
  - Added function plot_GIC

Thank you for reviewing our re-submission!

## Test environments
* local R installation (macOS 10.15.6), R 3.6.3 and 4.1.0
* local R installation (Ubuntu 20.04), R 4.0.5
* local R installation compiled without long double support (Ubuntu 20.04), R 4.0.3
* Ubuntu 16.04 (on travis-ci), R 4.0.2
* Windows Server x64 (on appveyor), R 4.0.3
* Win-builder (3.6.3, 4.0.3, and devel)

## R CMD check results (on local macOS: 3.6.3 and 4.1.0)

0 errors | 0 warnings | 1 note

* installed size is 12.3Mb
  sub-directories of 1Mb or more:
    * R:      2.6Mb
    * data:   9.0Mb
    
## Downstream dependencies

 There are currently no downstream dependencies for this package
