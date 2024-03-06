## Re-submission of MAPpoly package

This is a re-submission of the MAPpoly package. This version includes changes to address the issues that led to its archival on CRAN, specifically the failure in installing the package using `--use-LTO`. We have made corrections to ensure compatibility and successful installation with `--use-LTO`.

Changes in this submission include:

  - Addressed the `--use-LTO` installation failure issue

Thank you for reviewing our re-submission!

## Test environments
* local R installation (macOS Sonoma 14.3.1), R 4.3.3
   - Apple clang version 14.0.0 (clang-1400.0.29.202)
   - GNU Fortran (GCC) 12.2.0
* local R installation (Ubuntu 22.04), R 4.3.3
   - gcc (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
   - GNU Fortran (Ubuntu 9.4.0-1ubuntu1~20.04.2) 9.4.0
   - Compiled using the following configuration:
   
      ./configure \
       --prefix=/opt/R/${R_VERSION} \
       --enable-R-shlib \
       --enable-memory-profiling \
       --with-blas \
       --with-lapack\
       --enable-lto=R
       
    And installed using `R CMD INSTALL --use-LTO mappoly_0.4.1.tar.gz`
* local R installation (Ubuntu 22.04), R devel
   - Intel(R) oneAPI DPC++/C++ Compiler 2024.0.2 (2024.0.2.20231213)
   - ifx (IFX) 2024.0.2 20231213
   
* Win-builder (release and devel)


## R CMD check results 

0 errors | 0 warnings | 2 notes

 - on local macOS: R 4.3.3
   * installed size is 8.2Mb
     sub-directories of 1Mb or more:
       * R:      2.8Mb
       * data:   3.0Mb
   * GNU make is a SystemRequirements.

## Downstream dependencies

We checked 3 reverse dependencies (3 from CRAN + 0 from Bioconductor), comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages
