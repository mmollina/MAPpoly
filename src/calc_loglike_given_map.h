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
 File: calc_loglike_given_map.h
 
 Description: Set of functions to be used with software R
 
 Compute the log-likelihood using Hidden Markov Models (HMM) given a genetic map
 
 Function Written by Marcelo Mollinari.
 
 Contact: mmollina@ncsu.edu
 First version: Feb 04, 2020
 Last update:   Feb 04, 2020
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

RcppExport SEXP loglike_hmm(SEXP ploidyR,
                            SEXP genoR,
                            SEXP phPR,
                            SEXP phQR,
                            SEXP rfR,
                            SEXP verboseR);