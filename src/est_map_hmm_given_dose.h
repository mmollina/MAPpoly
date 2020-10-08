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

#ifndef _est_EST_H
#define _est_EST_H
#include <Rcpp.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include <math.h>
#include <R.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#include "combinatorial.h"

RcppExport SEXP  est_map_hmm(SEXP ploidyR, SEXP genoR, SEXP phPR, SEXP phQR,
                             SEXP rfR, SEXP verboseR, SEXP rf_limR, SEXP tol);

RcppExport SEXP est_map_hmm_highprec(SEXP ploidyR, SEXP genoR, SEXP phPR,
				     SEXP phQR, SEXP rfR, SEXP verboseR,
				     SEXP rf_limR, SEXP tolR);
#endif
