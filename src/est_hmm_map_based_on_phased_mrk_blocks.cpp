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
 File: est_hmm_map_based_on_phased_mrk_blocks.cpp
 
 Description: Set of functions to be used with software R
 
 Implements the methodology of Hidden Markov Models (HMM) to
 construct multipoint linkage maps in full-sib populations in
 autopolyploid species.
 
 Functions Written by Marcelo Mollinari.
 
 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2014
 Last update: Feb 18, 2016
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
#include "est_hmm_map_based_on_phased_mrk_blocks.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

RcppExport SEXP est_haplotype_map(SEXP ploidyR,
                                  SEXP n_marR,
                                  SEXP n_indR,
                                  SEXP haploR,
                                  SEXP emitR,
                                  SEXP rfR,
                                  SEXP verboseR,
                                  SEXP tolR,
                                  SEXP ret_H0R)
{
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mar = Rcpp::as<int>(n_marR);
  int n_ind = Rcpp::as<int>(n_indR);
  Rcpp::List haplo(haploR);
  Rcpp::List emit(emitR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  int ret_H0 = Rcpp::as<int>(ret_H0R);
  double tol = Rcpp::as<double>(tolR);
  
  //Initializing some variables
  int k, k1,  maxit = 1000, flag=0;
  double s, loglike=0.0, nr=0.0, temp=0.0;
  std::vector<double> rf_cur(rf.size());
  std::vector<double> term(n_ind);
  std::fill(term.begin(), term.end(), 0.0);
  
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mar << "\n" ;
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
    for(int i=0; i < n_mar; i++)
    {
      std::vector<double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
    }
  }
  
  //Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);
  
  //begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mar-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    //Initializing transition matrices
    std::vector< std::vector< std::vector<double> > > T;
    for(int i=0; i < n_mar-1; i++)
    {
      T.push_back(transition(m, rf_cur[i]));
    }
    //Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      {
        alpha[ind][0][j] = e[0][ind][j];
      }
      std::fill(beta[ind][n_mar-1].begin(), beta[ind][n_mar-1].end(), 1);
      //forward-backward
      for(k=1,k1=n_mar-2; k < n_mar; k++, k1--)
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
      if(ret_H0 == 0)
      {
        //Updating recombination fraction
        for(int k = 0; k < n_mar-1; k++)
        {
          vector<vector<double> > gamma(alpha[ind][k].size(), vector<double>(beta[ind][k+1].size()));
          s=0.0;
          int ngeni = alpha[ind][k].size();
          int ngenj = beta[ind][k+1].size();
          for(int i = 0; i < ngeni; i++)
          {
            for(int j = 0; j < ngenj; j++)
            {
              gamma[i][j] = alpha[ind][k][i] * beta[ind][k+1][j] *
                T[k][v[k][ind][i]][v[k+1][ind][j]] *
                T[k][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(i==0 && j==0) s = gamma[i][j];
              else s += gamma[i][j];
            }
          }
          for(int i=0; i < ngeni; i++)
          {
            for(int j=0; j < ngenj; j++)
            {
              nr=R[v[k][ind][i]][v[k+1][ind][j]] +
                R[v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(s > 0) // Verify theoretical implications of this condition
                rf[k] +=  nr * gamma[i][j]/s;
            }
          }
        }
      }
      //Termination
      for(int j=0; (unsigned)j < alpha[ind][n_mar-1].size(); j++)
      {
        term[ind] +=  alpha[ind][n_mar-1][j];
      }
    } // loop over individuals
    
    //Likelihood using a specific recombination fraction vector
    //Usually, this is used to compute LOD Score under H0: rf=0.5
    if(ret_H0 == 1)
    {
      //Loglike computation
      for(int i=0; (unsigned)i < alpha.size(); i++)
      {
        temp=0.0;
        for(int j=0; (unsigned)j < alpha[i][n_mar-1].size(); j++)
          temp += alpha[i][n_mar-1][j];
        if(temp > 0)
          loglike += log10(temp);
      }
      if(verbose)
        Rcpp::Rcout << "\n";
      List z = List::create(wrap(loglike), rf_cur);
      return(z);
    }
    
    // rescale
    for(int j=0; j<n_mar-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    // check convergence
    flag=0;
    for(int j=0; j < n_mar-1; j++)
    {
      if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
      {
        flag = 1;
        break;
      }
    }
    if(verbose)
    {
      Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
      for(int j = 0; j < n_mar-1; j++)
      {
        Rcpp::Rcout.precision(3);
        Rcpp::Rcout << std::fixed << rf[j] << " ";
      }
    }
    if(!flag) break;
  }//end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";
  
  //Loglike computation
  for(int i=0; (unsigned)i < alpha.size(); i++)
  {
    temp=0.0;
    for(int j=0; (unsigned)j < alpha[i][n_mar-1].size(); j++)
      temp += alpha[i][n_mar-1][j];
    if(temp > 0)
      loglike += log10(temp);
  }
  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}

RcppExport SEXP est_haplotype_map_highprec(SEXP ploidyR,
                                           SEXP n_marR,
                                           SEXP n_indR,
                                           SEXP haploR,
                                           SEXP emitR,
                                           SEXP rfR,
                                           SEXP verboseR,
                                           SEXP tolR,
                                           SEXP ret_H0R)
{
  //convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  int n_mar = Rcpp::as<int>(n_marR);
  int n_ind = Rcpp::as<int>(n_indR);
  Rcpp::List haplo(haploR);
  Rcpp::List emit(emitR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  int ret_H0 = Rcpp::as<int>(ret_H0R);
  double tol = Rcpp::as<double>(tolR);
  
  //Initializing some variables
  int k, k1,  maxit = 1000, flag=0;
  double s, loglike=0.0, nr=0.0;
  long double temp=0.0;
  std::vector<double> rf_cur(rf.size());
  std::vector<double> term(n_ind);
  std::fill(term.begin(), term.end(), 0.0);
  
  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mar << "\n" ;
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
  std::vector<std::vector<std::vector<long double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<long double> > > beta(n_ind);
  for(int ind=0; ind < n_ind; ind++)
  {
    for(int i=0; i < n_mar; i++)
    {
      std::vector<long double> temp3(v[i][ind].size()/2);
      alpha[ind].push_back(temp3);
      beta[ind].push_back(temp3);
    }
  }
  
  //Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);
  
  //begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mar-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    //Initializing transition matrices
    std::vector< std::vector< std::vector<double> > > T;
    for(int i=0; i < n_mar-1; i++)
    {
      T.push_back(transition(m, rf_cur[i]));
    }
    //Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      for(int j=0; (unsigned)j < e[0][ind].size(); j++)
      {
        alpha[ind][0][j] = e[0][ind][j];
      }
      std::fill(beta[ind][n_mar-1].begin(), beta[ind][n_mar-1].end(), 1);
      //forward-backward
      for(k=1,k1=n_mar-2; k < n_mar; k++, k1--)
      {
        std::vector<long double> temp4 (v[k][ind].size()/2);
        temp4 = forward_emit_highprec(m, alpha[ind][k-1], v[k-1][ind], v[k][ind], e[k][ind], T[k-1]);
        for(int j=0; (unsigned)j < temp4.size(); j++)
        {
          alpha[ind][k][j]=temp4[j];
        }
        std::vector<long double> temp5 (v[k1][ind].size()/2);
        temp5=backward_emit_highprec(m, beta[ind][k1+1], v[k1][ind], v[k1+1][ind], e[k1+1][ind], T[k1]);
        for(int j=0; (unsigned)j < temp5.size(); j++)
        {
          beta[ind][k1][j]=temp5[j];
        }
      }
      if(ret_H0 == 0)
      {
        //Updating recombination fraction
        for(int k = 0; k < n_mar-1; k++)
        {
          vector<vector<long double> > gamma(alpha[ind][k].size(), vector<long double>(beta[ind][k+1].size()));
          s=0.0;
          int ngeni = alpha[ind][k].size();
          int ngenj = beta[ind][k+1].size();
          for(int i = 0; i < ngeni; i++)
          {
            for(int j = 0; j < ngenj; j++)
            {
              gamma[i][j] = alpha[ind][k][i] * beta[ind][k+1][j] *
                T[k][v[k][ind][i]][v[k+1][ind][j]] *
                T[k][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(i==0 && j==0) s = gamma[i][j];
              else s += gamma[i][j];
            }
          }
          for(int i=0; i < ngeni; i++)
          {
            for(int j=0; j < ngenj; j++)
            {
              nr=R[v[k][ind][i]][v[k+1][ind][j]] +
                R[v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(s > 0) // Verify theoretical implications of this condition
                rf[k] +=  nr * gamma[i][j]/s;
            }
          }
        }
      }
      //Termination
      for(int j=0; (unsigned)j < alpha[ind][n_mar-1].size(); j++)
      {
        term[ind] +=  alpha[ind][n_mar-1][j];
      }
    } // loop over individuals
    
    //Likelihood using a specific recombination fraction vector
    //Usually, this is used to compute LOD Score under H0: rf=0.5
    if(ret_H0 == 1)
    {
      //Loglike computation
      for(int i=0; (unsigned)i < alpha.size(); i++)
      {
        temp=0.0;
        for(int j=0; (unsigned)j < alpha[i][n_mar-1].size(); j++)
          temp += alpha[i][n_mar-1][j];
        if(temp > 0)
          loglike += log10(temp);
      }
      if(verbose)
        Rcpp::Rcout << "\n";
      List z = List::create(wrap(loglike), rf_cur);
      return(z);
    }
    
    // rescale
    for(int j=0; j<n_mar-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    // check convergence
    flag=0;
    for(int j=0; j < n_mar-1; j++)
    {
      if(fabs(rf[j] - rf_cur[j]) > tol*(rf_cur[j]+tol*100.0))
      {
        flag = 1;
        break;
      }
    }
    if(verbose)
    {
      Rcpp::Rcout << "\t\n Iter: " << it+1 << "\t";
      for(int j = 0; j < n_mar-1; j++)
      {
        Rcpp::Rcout.precision(3);
        Rcpp::Rcout << std::fixed << rf[j] << " ";
      }
    }
    if(!flag) break;
  }//end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";
  
  //Loglike computation
  for(int i=0; (unsigned)i < alpha.size(); i++)
  {
    temp=0.0;
    for(int j=0; (unsigned)j < alpha[i][n_mar-1].size(); j++)
      temp += alpha[i][n_mar-1][j];
    if(temp > 0)
      loglike += log10(temp);
  }
  if(verbose)
    Rcpp::Rcout << "\n";
  List z = List::create(wrap(loglike), rf);
  return(z);
}


