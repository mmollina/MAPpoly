/*
 MAPpoly: a package to construct genetic maps in autopolyploids
Copyright (C) 2014-2018 Marcelo Mollinari

This file is part of MAPpoly.

MAPpoly is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
File: pairwise_estimation.h

Description: Set of functions to be used with software R

Implements the methodology of Hidden Markov Models (HMM) to
construct multipoint linkage maps in full-sib populations in
autopolyploid species

Functions Written by Marcelo Mollinari.

Contact: mmollina@ncsu.edu
First version:       2014
Last update: Feb 18, 2016
*/

#ifndef PAIRWISE_ESTIMATION_H
#define PAIRWISE_ESTIMATION_H
#include <Rcpp.h>

double twopt_likelihood_dosage(double rf, int m, int n_ind,
			       int dP, int dQ,
			       Rcpp::NumericVector dk,
			       Rcpp::NumericVector dk1,
			       Rcpp::NumericVector gen_1,
			       Rcpp::NumericVector gen_2,
			       Rcpp::NumericMatrix count_mat);

RcppExport SEXP pairwise_rf_estimation(SEXP m_R,
				       SEXP mrk_pairs_R,
				       SEXP geno_R,
				       SEXP dP_R,
				       SEXP dQ_R,
				       SEXP count_cache_R,
				       SEXP tol_R);
#endif