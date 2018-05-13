#include <Rcpp.h>
#include "genotypic_counts.h"
#include "combinatorial.h"

using namespace Rcpp;

RcppExport SEXP get_counts_one_parent_cpp(SEXP ploidyR, SEXP gen_par_mk1R,
					  SEXP gen_par_mk2R, SEXP gen_prog_mk1R,
					  SEXP gen_prog_mk2R, SEXP countsR)
{

  int ploidy = Rcpp::as<int>(ploidyR);
  int gen_prog_mk1 = Rcpp::as<int>(gen_prog_mk1R);
  int gen_prog_mk2 = Rcpp::as<int>(gen_prog_mk2R);
  std::vector<int> gen_par_mk1 = Rcpp::as<std::vector<int> >(gen_par_mk1R);
  std::vector<int> gen_par_mk2 = Rcpp::as<std::vector<int> >(gen_par_mk2R);
  std::vector<int> counts = Rcpp::as<std::vector<int> >(countsR);

  counts=boolean_lexicographic_k_choose_m_and_collapse(ploidy,
						       gen_par_mk1,
						       gen_par_mk2,
						       gen_prog_mk1,
						       gen_prog_mk2);
  List z  = List::create(ploidy,
			 gen_par_mk1,
			 gen_par_mk2,
			 gen_prog_mk1,
			 gen_prog_mk2,
			 counts) ;
  return z ;
}

