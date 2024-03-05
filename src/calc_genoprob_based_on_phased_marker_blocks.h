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

#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "combinatorial.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

List calc_genprob_haplo_cpp(int m, 
                            int n_mrk, 
                            int n_ind, 
                            List haplo, 
                            List emit, 
                            NumericVector rf, 
                            std::vector<long double> probs, 
                            int verbose);

List calc_genprob_haplo_highprec(int m, 
                                 int n_mrk, 
                                 int n_ind, 
                                 List haplo, 
                                 List emit, 
                                 NumericVector rf, 
                                 std::vector<long double> probs, 
                                 int verbose);