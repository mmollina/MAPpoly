/*
 Polymap: a package to construct genetic maps in autopolyploids
Copyright (C) 2014-2016 Marcelo Mollinari

This file is part of Polymap.

Polymap is free software: you can redistribute it and/or modify
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
File: est_map_hmm.cpp

Description: Set of functions to be used with software R

Implements the methodology of Hidden Markov Models (HMM) to
construct multipoint linkage maps in full-sib populations in
autopolyploid species

Functions Written by Marcelo Mollinari.

Contact: mmollina@ncsu.edu
First version:       Dec 6, 2019
Last update: Dec 6, 2019
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

RcppExport SEXP calc_genprob_haplo(SEXP ploidyR,
                                   SEXP n_mrkR,
                                   SEXP n_indR,
                                   SEXP haploR,
                                   SEXP emitR,
                                   SEXP rfR,
                                   SEXP probsR,
                                   SEXP verboseR);