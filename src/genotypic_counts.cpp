/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2020 Marcelo Mollinari
 
 This file is part of MAPpoly.
 
 MAPpoly is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 For a copy of the GNU General Public License, please visit
 <http://www.gnu.org/licenses/>.
 */
/*
 Functions Written by Marcelo Mollinari.
 
 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 */
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

