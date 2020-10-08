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
 File: calc_loglike_given_map.cpp
 
 Description: Set of functions to be used with software R
 
 Compute the log-likelihood using Hidden Markov Models (HMM) given a genetic map
 
 Function Written by Marcelo Mollinari.
 
 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
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
#include "calc_loglike_given_map.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

RcppExport SEXP loglike_hmm(SEXP ploidyR,
                            SEXP genoR,
                            SEXP phPR,
                            SEXP phQR,
                            SEXP rfR,
                            SEXP verboseR)
{
  //*convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  Rcpp::NumericMatrix geno(genoR); //(n.ind x n.col)
  Rcpp::List ph1(phPR);
  Rcpp::List ph2(phQR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  
  //*Initializing some variables
  int g = nChoosek(m, m/2);
  int n_mar = geno.ncol(); // markers are disposed in columns
  int n_ind = geno.nrow(); // individuals are disposed in rows
  double loglike=0.0, temp=0.0;
  
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mar << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
    for(int j = 0; j < n_mar-1; j++)
    {
      Rcpp::Rcout.precision(3);
      Rcpp::Rcout << std::fixed << rf[j] << " ";
    }
  }
  //*Initializing v: states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  for(int i=0; i < n_mar; i++)
  {
    std::vector<int> temp1 = as<vector<int> > (ph1[i]);
    std::vector<int> temp2 = as<vector<int> > (ph2[i]);
    v.push_back(index_func(m, temp1, temp2));
  }
  //*Initializing alpha 
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  for(int k=0; k < n_ind; k++)
  {
    for(int i=0; i < n_mar; i++)
    {
      std::vector<double> temp3(v[i][geno(k,i)].size()/2);
      alpha[k].push_back(temp3);
    }
  }
  
  //*Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);
  
  //*Initializing transition matrices
  std::vector< std::vector< std::vector<double> > > T;
  for(int i=0; i < n_mar-1; i++)
  {
    T.push_back(transition(m, rf[i]));
  }
  //*Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    std::fill(alpha[ind][0].begin(), alpha[ind][0].end(), 1.0/(g*g));
    //*forward
    for(int k=1; k < n_mar; k++)
    {
      std::vector<double> temp4 (v[k][geno(ind,k)].size()/2);
      temp4 = forward(m, alpha[ind][k-1], v[k-1][geno(ind,k-1)], v[k][geno(ind,k)], T[k-1]);
      for(int j=0; (unsigned)j < temp4.size(); j++)
      {
        alpha[ind][k][j]=temp4[j];
      }
    }
  } //* loop over individuals
  for(int i=0; (unsigned)i < alpha.size(); i++)
  {
    temp=0.0;
    for(int j=0; (unsigned)j < alpha[i][n_mar-1].size(); j++)
    {
      temp += alpha[i][n_mar-1][j];
    }
    if(temp > 0)
      loglike += log10(temp);
  }
  if(verbose) Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}
