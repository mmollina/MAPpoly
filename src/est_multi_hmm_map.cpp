/*
 MAPpoly: a package to construct genetic maps in autopolyploids
 Copyright (C) 2014-2022 Marcelo Mollinari

 This file is part of MAPpoly.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 */

/*
 File: est_hmm_map.cpp

 Description:

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version:       2022
 Last update: Oct 06, 2022
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

RcppExport SEXP est_hmm_map(SEXP ploidy1_idR,
                            SEXP ploidy2_idR,
                            SEXP max_ploidy_idR,
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
  Rcpp::NumericVector p1(ploidy1_idR);
  Rcpp::NumericVector p2(ploidy2_idR);
  int max_ploidy_id = Rcpp::as<int>(max_ploidy_idR);
  int n_mar = Rcpp::as<int>(n_marR);
  int n_ind = Rcpp::as<int>(n_indR);
  Rcpp::List haplo(haploR);
  Rcpp::List emit(emitR);
  Rcpp::NumericVector rf(rfR);
  int verbose = Rcpp::as<int>(verboseR);
  int ret_H0 = Rcpp::as<int>(ret_H0R);
  double tol = Rcpp::as<double>(tolR);
  std::vector<int> pl{2,4,6};
  //Initializing some variables
  int k, k1,  maxit = 1000, flag=0;
  double s, loglike=0.0, nr=0.0, temp=0.0;
  std::vector<double> rf_cur(rf.size());
  std::vector<double> term(n_ind);
  std::fill(term.begin(), term.end(), 0.0);

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
  std::vector< std::vector< std::vector<double> > > R;
  for(int j=0; j <= max_ploidy_id; j++)
    R.push_back(rec_num(pl[j]));

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
    std::vector< std::vector< std::vector< std::vector<double> > > >T;
    for(int j=0; j <= max_ploidy_id; j++)
    {
      std::vector< std::vector< std::vector<double> > > Ttemp;
      for(int i=0; i < n_mar-1; i++)
      {
        Ttemp.push_back(transition(pl[j], rf_cur[i]));
      }
      T.push_back(Ttemp);
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
        temp4 = forward_emit(alpha[ind][k-1], v[k-1][ind], v[k][ind], e[k][ind], T[p1[ind]][k-1], T[p2[ind]][k-1]);
        for(int j=0; (unsigned)j < temp4.size(); j++)
        {
          alpha[ind][k][j]=temp4[j];
        }
        std::vector<double> temp5 (v[k1][ind].size()/2);
        temp5=backward_emit(beta[ind][k1+1], v[k1][ind], v[k1+1][ind], e[k1+1][ind], T[p1[ind]][k1], T[p2[ind]][k1]);
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
                T[p1[ind]][k][v[k][ind][i]][v[k+1][ind][j]] *
                T[p2[ind]][k][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
              if(i==0 && j==0) s = gamma[i][j];
              else s += gamma[i][j];
            }
          }
          for(int i=0; i < ngeni; i++)
          {
            for(int j=0; j < ngenj; j++)
            {
              nr=R[p1[ind]][v[k][ind][i]][v[k+1][ind][j]] +
                R[p2[ind]][v[k][ind][i+ngeni]][v[k+1][ind][j+ngenj]];
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
