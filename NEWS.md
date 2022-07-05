# MAPpoly version 0.3.1
 - Minor modifications to the test functions
 - Changed 'smacof' package to suggested  

# MAPpoly version 0.3.0
 - Added function est_pairwise_rf2 to avoid memory overflow in  personal computers when estimating recombination fraction in large number of markers.
 - Added Vignette
 - Added function 'export_qtlpoly'
 - Added function 'filter_individuals'
 - Changed the name of function 'calc_homoprob' to 'calc_homologprob'
 
# MAPpoly version 0.2.3
 - Added function plot_GIC
 - Updated functions est_rf_hmm_sequential and read_fitpoly
 - Verification for high-prec system support
 - Added LazyDataCompression: xz in DESCRIPTION file
 
# MAPpoly version 0.2.2
  - Suggested packages are used conditionally, following §1.1.3.1 of 'Writing R Extensions'.

# MAPpoly version 0.2.1
  - Added estimation of two-point recombination fraction using genotype probabilities
  - Added 'read_fitpoly' function
  - New graphical feature in function 'rf_snp_filt'
  - Added option to plot large recombination fraction matrix using cell aggregation for a given order
  - Added estimation of map distance using the MDSMap package procedure (projects the MDS order onto a single dimension using principal curves)
  
  
