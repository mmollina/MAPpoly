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
 File: est_map_hmm_given_dose.cpp

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
#include "est_map_hmm_given_dose.h"
#include "hmm_elements.h"
using namespace std;
using namespace Rcpp;

RcppExport SEXP est_map_hmm(SEXP ploidyR,
                            SEXP genoR,
                            SEXP phPR,
                            SEXP phQR,
                            SEXP rfR,
                            SEXP verboseR,
                            SEXP rf_limR,
                            SEXP tolR)
{
  //*convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  Rcpp::NumericMatrix geno(genoR); //(n.ind x n.col)
  Rcpp::List ph1(phPR);
  Rcpp::List ph2(phQR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  double rf_lim = Rcpp::as<double>(rf_limR);
  double tol = Rcpp::as<double>(tolR);

  //*Initializing some variables
  int g = nChoosek(m, m/2), k, k1,  maxit = 1000, flag=0;
  int n_mar = geno.ncol(); // markers are disposed in columns
  int n_ind = geno.nrow(); // individuals are disposed in rows
  double s, loglike=0.0, nr=0.0, temp=0.0;
  std::vector<double> rf_cur(rf.size());

  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tRec. Frac. Limit: " << rf_lim << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mar << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
    Rcpp::Rcout << "\t\n Init. values:\t";
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
  //*Initializing alpha and beta
  std::vector<std::vector<std::vector<double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<double> > > beta(n_ind);
  for(int k=0; k < n_ind; k++)
  {
    for(int i=0; i < n_mar; i++)
    {
      std::vector<double> temp3(v[i][geno(k,i)].size()/2);
      alpha[k].push_back(temp3);
      beta[k].push_back(temp3);
    }
  }

  //*Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);

  //*begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //*Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mar-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    //*Initializing transition matrices
    std::vector< std::vector< std::vector<double> > > T;
    for(int i=0; i < n_mar-1; i++)
    {
      T.push_back(transition(m, rf_cur[i]));
    }
    //*Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      std::fill(alpha[ind][0].begin(), alpha[ind][0].end(), 1.0/(g*g));
      std::fill(beta[ind][n_mar-1].begin(), beta[ind][n_mar-1].end(), 1);
      //*forward-backward
      for(k=1,k1=n_mar-2; k < n_mar; k++, k1--)
      {
        std::vector<double> temp4 (v[k][geno(ind,k)].size()/2);
        temp4 = forward(m, alpha[ind][k-1], v[k-1][geno(ind,k-1)], v[k][geno(ind,k)], T[k-1]);
        for(int j=0; (unsigned)j < temp4.size(); j++)
        {
          alpha[ind][k][j]=temp4[j];
        }
        std::vector<double> temp5 (v[k1][geno(ind,k1)].size()/2);
        temp5=backward(m, beta[ind][k1+1], v[k1][geno(ind,k1)], v[k1+1][geno(ind,k1+1)], T[k1]);
        for(int j=0; (unsigned)j < temp5.size(); j++)
        {
          beta[ind][k1][j]=temp5[j];
        }
      }
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
              T[k] [v[k][geno(ind,k)][i]] [v[k+1][geno(ind,k+1)][j]] *
              T[k] [v[k][geno(ind,k)][i+ngeni]] [v[k+1][geno(ind,k+1)][j+ngenj]];
            if(i==0 && j==0) s = gamma[i][j];
            else s += gamma[i][j];
          }
        }
        for(int i=0; i < ngeni; i++)
        {
          for(int j=0; j < ngenj; j++)
          {
            nr=R[v[k][geno(ind,k)][i]][v[k+1][geno(ind,k+1)][j]] +
              R[v[k][geno(ind,k)][i+ngeni]][v[k+1][geno(ind,k+1)][j+ngenj]];
            if(s > 0)
              rf[k] +=  (double)nr * gamma[i][j]/s;
          }
        }
      }
    } //* loop over individuals
    if(tol>=0.99)
    {
      for(int j=0; j < n_mar-1; j++)
      {
        rf[j] = rf_cur[j];
      }
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
    //* rescale
    for(int j=0; j<n_mar-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    flag=0;
    //* check convergence
    for(int j=0; j < n_mar-1; j++)
    {
      if(rf[j] > rf_lim)
      {
        flag = 0;
        break; // if any element in rf_cur is greater than rf_lim, return the current values
      }
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
  }//*end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";
  //*Loglike computation

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


/*
Estimates a genetic map usinh HMM technology given
an order and a linkage phase configuration.
High precision version
*/
RcppExport SEXP est_map_hmm_highprec(SEXP ploidyR,
				     SEXP genoR,
				     SEXP phPR,
				     SEXP phQR,
				     SEXP rfR,
				     SEXP verboseR,
				     SEXP rf_limR,
				     SEXP tolR)
{
  //*convert input to C++ types
  int m = Rcpp::as<int>(ploidyR);
  Rcpp::NumericMatrix geno(genoR); //(n.ind x n.col)
  Rcpp::List ph1(phPR);
  Rcpp::List ph2(phQR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  double rf_lim = Rcpp::as<double>(rf_limR);
  double tol = Rcpp::as<double>(tolR);

  //*Initializing some variables
  int g = nChoosek(m, m/2), k, k1,  maxit = 1000, flag=0;
  int n_mar = geno.ncol(); // markers are disposed in columns
  int n_ind = geno.nrow(); // individuals are disposed in rows
  double nr=0.0;
  long double s, loglike=0.0, temp=0.0;
  std::vector<double> rf_cur(rf.size());

  if(verbose)
  {
    Rcpp::Rcout << "\tPloidy level: " << m << "\n" ;
    Rcpp::Rcout << "\tRec. Frac. Limit: " << rf_lim << "\n" ;
    Rcpp::Rcout << "\tNumber of markers: " << n_mar << "\n" ;
    Rcpp::Rcout << "\tNumber of individuals: " << n_ind << "\n" ;
    Rcpp::Rcout << "\t\n Init. values:\t";
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
  //*Initializing alpha and beta
  std::vector<std::vector<std::vector<long double> > > alpha(n_ind);
  std::vector<std::vector<std::vector<long double> > > beta(n_ind);
  for(int k=0; k < n_ind; k++)
  {
    for(int i=0; i < n_mar; i++)
    {
      std::vector<long double> temp3(v[i][geno(k,i)].size()/2);
      alpha[k].push_back(temp3);
      beta[k].push_back(temp3);
    }
  }

  //*Initializing recombination number matrix
  std::vector< std::vector<double> > R;
  R = rec_num(m);

  //*begin EM algorithm
  for(int it=0; it<maxit; it++)
  {
    //*Initializing recombination fraction vector for Baum-Welch
    for(int j=0; j<n_mar-1; j++)
    {
      rf_cur[j] = rf[j];
      rf[j] = 0.0;
    }
    //*Initializing transition matrices
    std::vector< std::vector< std::vector<double> > > T;
    for(int i=0; i < n_mar-1; i++)
    {
      T.push_back(transition(m, rf_cur[i]));
    }
    //*Loop over all individuals
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      std::fill(alpha[ind][0].begin(), alpha[ind][0].end(), 1.0/(g*g));
      std::fill(beta[ind][n_mar-1].begin(), beta[ind][n_mar-1].end(), 1);
      //*forward-backward
      for(k=1,k1=n_mar-2; k < n_mar; k++, k1--)
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
              T[k] [v[k][geno(ind,k)][i]] [v[k+1][geno(ind,k+1)][j]] *
              T[k] [v[k][geno(ind,k)][i+ngeni]] [v[k+1][geno(ind,k+1)][j+ngenj]];
            if(i==0 && j==0) s = gamma[i][j];
            else s += gamma[i][j];
          }
        }
        for(int i=0; i < ngeni; i++)
        {
          for(int j=0; j < ngenj; j++)
          {
            nr=R[v[k][geno(ind,k)][i]][v[k+1][geno(ind,k+1)][j]] +
              R[v[k][geno(ind,k)][i+ngeni]][v[k+1][geno(ind,k+1)][j+ngenj]];
            if(s > 0)
              rf[k] +=  (long double)nr * gamma[i][j]/s;
          }
        }
      }
    } //* loop over individuals
    if(tol>=0.99)
    {
      for(int j=0; j < n_mar-1; j++)
      {
        rf[j] = rf_cur[j];
      }
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
    //* rescale
    for(int j=0; j<n_mar-1; j++)
    {
      rf[j] /= (double)n_ind;
      if(rf[j] < tol/100.0) rf[j] = tol/100.0;
      else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
    }
    flag=0;
    //* check convergence
    for(int j=0; j < n_mar-1; j++)
    {
      if(rf[j] > rf_lim)
      {
        flag = 0;
        break; // if any element in rf_cur is greater than rf_lim, return the current values
      }
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
  }//*end of EM algorithm
  if(flag && verbose) Rcpp::Rcout << "Didn't converge!\n";
  //*Loglike computation

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

