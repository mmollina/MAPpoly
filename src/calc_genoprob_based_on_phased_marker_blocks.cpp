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
 File: calc_genoprob_based_on_phased_marker_blocks.cpp
 
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
 Last update: 6 Dec,  2019
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
                                   SEXP verboseR)
{
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mrk = Rcpp::as<int>(n_mrkR);
  int n_ind = Rcpp::as<int>(n_indR);
  Rcpp::List haplo(haploR);
  Rcpp::List emit(emitR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  std::vector<long double> probs = Rcpp::as<std::vector<long double> >(probsR);
  
  //Initializing some variables
  int g = nChoosek(m, m/2), k, k1, count = 0;
  std::vector<double> term(n_ind);
  std::fill(term.begin(), term.end(), 0.0);
  
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mrk << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
  }
  //Initializing v: states hmm should visit for each marker
  //Initializing e: emission probabilities associated to the states hmm should visit for each marker
  std::vector<std::vector<std::vector<int> > > v;
  std::vector<std::vector<std::vector<double> > > e;
  for(int i=0; i < haplo.size(); i++) // i: number of markers
  {
    Rcpp::List haplo_temp(haplo(i)); //states hmm should visit for marker i
    Rcpp::List emit_temp(emit(i)); //emission probs. for states hmm should visit for marker i
    std::vector<std::vector<int> > v1;
    std::vector<std::vector<double> > e1;
    for(int j=0; j < haplo_temp.size(); j++) //iterate for all j individuals
    {
      Rcpp::NumericMatrix M_temp = haplo_temp(j);
      Rcpp::NumericVector E_temp = emit_temp(j);
      std::vector<int> v2 = Rcpp::as<std::vector<int> >(M_temp);
      std::vector<double> e2 = Rcpp::as<std::vector<double> >(E_temp);
      v1.push_back(v2);
      e1.push_back(e2);
    }
    v.push_back(v1);
    e.push_back(e1);
  }
  //Initializing alpha and beta
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<double> > > beta(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mrk; i++)
    {
      std::vector<double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
    }
  }
  //Initializing transition matrices
  std::vector< std::vector< std::vector<double> > > T;
  for(int i=0; i < n_mrk-1; i++)
  {
    T.push_back(transition(m, rf[i]));
  }
  //Loop over all individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    R_CheckUserInterrupt();
    //std::fill(alpha[ind][0].begin(), alpha[ind][0].end(), 1.0/(g*g));
    for(int j=0; (unsigned)j < e[0][ind].size(); j++)
    {
      alpha[ind][0][j] = e[0][ind][j];
    }
    std::fill(beta[ind][n_mrk-1].begin(), beta[ind][n_mrk-1].end(), 1);
    //forward-backward
    for(k=1,k1=n_mrk-2; k < n_mrk; k++, k1--)
    {
      std::vector<double> temp4 (v[k][ind].size()/2);
      temp4 = forward_emit(m, alpha[ind][k-1], v[k-1][ind], v[k][ind], e[k][ind], T[k-1]);
      for(int j=0; (unsigned)j < temp4.size(); j++)
      {
        alpha[ind][k][j]=temp4[j];
      }
      std::vector<double> temp5 (v[k1][ind].size()/2);
      temp5=backward_emit(m, beta[ind][k1+1], v[k1][ind], v[k1+1][ind], e[k1+1][ind], T[k1]);
      for(int j=0; (unsigned)j < temp5.size(); j++)
      {
        beta[ind][k1][j]=temp5[j];
      }
    }
    //Termination
    for(int j=0; (unsigned)j < alpha[ind][n_mrk-1].size(); j++)
    {
      term[ind] +=  alpha[ind][n_mrk-1][j];
    }
  } // loop over individuals
  for(int ind=0; ind < n_ind; ind++)
  {
    std::vector<double> gamma(g*g);
    for(int k=0; k < n_mrk; k++)
    {
      std::fill(gamma.begin(), gamma.end(), 0);
      int s = v[k][ind].size();
      long double w = 0.0;
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        w += alpha[ind][k][j]*beta[ind][k][j];
      for(int j=0; (unsigned)j < alpha[ind][k].size(); j++)
        gamma[v[k][ind][j]*g + v[k][ind][j+s/2]] = alpha[ind][k][j]*beta[ind][k][j]/w;
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
