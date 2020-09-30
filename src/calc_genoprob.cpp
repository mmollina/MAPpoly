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
 File: calc_genoprob.cpp

 Description: Set of functions to be used with software R

 Implements the methodology of Hidden Markov Models (HMM) to
 construct multipoint linkage maps in full-sib populations in
 autopolyploid species

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2014
 Last update: Oct 16, 2017
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
#include "est_map_hmm_given_dose.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

RcppExport SEXP calc_genoprob(SEXP ploidyR,
                              SEXP genoR,
                              SEXP phPR,
                              SEXP phQR,
                              SEXP rfR,
                              SEXP probsR,
                              SEXP verboseR)
{
  //*convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  Rcpp::NumericMatrix geno(genoR); //(n.ind x n.col)
  Rcpp::List ph1(phPR);
  Rcpp::List ph2(phQR);
  Rcpp::NumericVector rf(rfR);
  std::vector<long double> probs = Rcpp::as<std::vector<long double> >(probsR);
  int verbose = Rcpp::as<int>(verboseR);

  //*Initializing some variables
  int g = nChoosek(m, m/2), k, k1;
  int n_mrk = geno.ncol(); // markers are disposed in columns
  int n_ind = geno.nrow(); // individuals are disposed in rows
  int count = 0;
  std::vector<double> rf_cur(rf.size());
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n\t" ;
  }
  //*Initializing v: states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  for(int i=0; i < n_mrk; i++)
  {
    std::vector<int> temp1 = as<vector<int> > (ph1[i]);
    std::vector<int> temp2 = as<vector<int> > (ph2[i]);
    v.push_back(index_func(m, temp1, temp2));
  }
  //*Initializing alpha and beta
  std::vector<std::vector<std::vector<long double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<long double> > > beta(n_ind);
  for(int k=0; k < n_ind; k++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<long double> temp3(v[i][geno(k,i)].size()/2);
      alpha[k].push_back(temp3);
      beta[k].push_back(temp3);
    }
  }

  //*Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);
  //*Initializing transition matrices
  std::vector< std::vector< std::vector<double> > > T;
  for(int i=0; i < n_mrk-1; i++)
  {
    T.push_back(transition(m, rf[i]));
  }
  //*Loop over all individuals

  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    if(verbose){
      Rcpp::Rcout << ".";
      if((ind+1)%50 == 0) Rcpp::Rcout << "\n\t";
    }
    std::fill(alpha[ind][0].begin(), alpha[ind][0].end(), 1.0/(g*g));
    std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1);
    //*forward-backward
    for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
    {
      std::vector<long double> temp4 (v[k][geno(ind,k)].size()/2);
      temp4 = forward_highprec(m, alpha[ind][k-1], v[k-1][geno(ind,k-1)], v[k][geno(ind,k)], T[k-1]);
      for(int j=0; (unsigned)j < temp4.size(); j++)
      {
        alpha[ind][k][j]=temp4[j];
      }
      std::vector<long double> temp5 (v[k1][geno(ind,k1)].size()/2);
      temp5=backward_highprec(m, beta[ind][k1+1], v[k1][geno(ind,k1)], v[k1+1][geno(ind,k1+1)], T[k1]);
      for(int j=0; (unsigned)j < temp5.size(); j++)
      {
        beta[ind][k1][j]=temp5[j];
      }
    }
  } //* loop over individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    std::vector<double> gamma(g*g);
    for(int k=0; k < n_mrk; k++)
    {
      std::fill(gamma.begin(), gamma.end(), 0);
      int s = v[k][geno(ind,k)].size();
      long double w = 0.0;
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        w += alpha[ind][k][j]*beta[ind][k][j];
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        gamma[v[k][geno(ind,k)][j]*g + v[k][geno(ind,k)][j+s/2]] = alpha[ind][k][j]*beta[ind][k][j]/w;
      for(int j=0; j < g*g; j++)
      {
        probs[count] = gamma[j];
        count++;
      }
    }
  }
  List z  = List::create(probs);
  return z ;
}

