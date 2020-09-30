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
 File: est_map_hmm_given_prior.cpp

 Description: Set of functions to be used with software R

 Implements the methodology of Hidden Markov Models (HMM) to
 construct multipoint linkage maps in full-sib populations in
 autopolyploid species

 Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version: Dec 19, 2013
 Last update: Sep 24, 2020
 */

#include <Rcpp.h>
#include "hmm_elements.h"
#include "combinatorial.h"
#include <math.h>
#include <algorithm>
#include <new>
#define TOL 0
#define THRESHOLD 0.01

using namespace std;
using namespace Rcpp;

const int MAX_SIZE = 853776;
int pre_calc_n_rec_1[MAX_SIZE] = {};
int pre_calc_n_rec_2[MAX_SIZE] = {};

void setup_pre_calc_n_rec_cache(int m, int gam)
{
  for(int i = 0; i < gam; ++i)
  {
    for(int j = 0; j < gam; ++j)
    {
      pre_calc_n_rec_2[(i*gam)+j] = pre_calc_n_rec_1[(i*gam)+j] = n_rec_given_genk_and_k1(m, i+1, j+1);
    }
  }
}

/* Function: alpha_ai_dot

 pre_calc_n_rec_1           pre_calc_n_rec_2
          --         --              --         --
 l1 = 0  |0 1 1 1 1 2|     l2 = 0   |0 1 1 1 1 2|
 l1 = 1  |1 0 1 1 2 1|     l2 = 1   |1 0 1 1 2 1|
 l1 = 2  |1 1 0 2 1 1|     l2 = 2   |1 1 0 2 1 1|
 l1 = 3  |1 1 2 0 1 1|     l2 = 3   |1 1 2 0 1 1|
 l1 = 4  |1 2 1 1 0 1|     l2 = 4   |1 2 1 1 0 1|
 l1 = 5  |2 1 1 1 1 0|     l2 = 5   |2 1 1 1 1 0|
         --         --              --         --

 Given the number of recombination provided by the two
 matrices above (tetraploid example), function
 pre_calc_prob returns Pr(pk+1|pk). 
 
 For example, l1=0 and l2=0 --> function pre_calc_prob returns P(pk+1 = {*}|pk={1,2,3,4}),
 i.e.
           --                      --
 1 2 3 4   |x x x x x x ....x x x x |  <-- l1=0, l2=0;
 1 2 3 5   |...                  ...|
            .
            .
            .
 */

std::vector <double> alpha_ai_dot(std::vector<double>& pre_calc_prob, int m, int gam, int gam_pow_2, int l1, int l2)
{
  std::vector<double> a(gam_pow_2);
  a[0]=pre_calc_prob[(1+m/2)*pre_calc_n_rec_1[0+(l1*gam)]+pre_calc_n_rec_2[0+(l2*gam)]];
  for(int j=1; j < gam; j++)
  {
    if(pre_calc_n_rec_2[j+(l2*gam)]!=pre_calc_n_rec_2[j-1+(l2*gam)])
    {
      a[j]=pre_calc_prob[(1+m/2)*pre_calc_n_rec_1[0+(l1*gam)]+pre_calc_n_rec_2[j+(l2*gam)]];
    }
    else
    {
      a[j]=a[j-1];
    }
  }
  for(int i=1; i < gam; i++)
  {
    if(pre_calc_n_rec_1[i+(l1*gam)]!=pre_calc_n_rec_1[i-1+(l1*gam)])
    {
      a[i*gam]=pre_calc_prob[(1+m/2)*pre_calc_n_rec_1[i+(l1*gam)]+pre_calc_n_rec_2[0+(l2*gam)]];
      for(int j=1; j < gam; j++)
      {
        if(pre_calc_n_rec_2[j+(l2*gam)]!=pre_calc_n_rec_2[j-1+(l2*gam)])
        {
          a[(i*gam)+j]=pre_calc_prob[(1+m/2)*pre_calc_n_rec_1[i+(l1*gam)]+pre_calc_n_rec_2[j+(l2*gam)]];
        }
        else
        {
          a[(i*gam)+j]=a[(i*gam)+j-1];
        }
      }
    }
    else
    {
      for(int j1 = (i-1)*gam; j1 < i*gam; j1++)
      {
        a[j1+gam] = a[j1];
      }
    }
  }
  return a;
}

std::vector <double> nrec(int m, int gam, int gam_pow_2, int l1, int l2)
{
  std::vector<double> rec(gam_pow_2);
  for(int i=0; i < gam; i++)
  {
    for(int j=0; j < gam; j++)
    {
      rec[(i*gam)+j]=(pre_calc_n_rec_1[i+(l1*gam)]+pre_calc_n_rec_2[j+(l2*gam)])/(double)m;
    }
  }
  return rec;
}

RcppExport SEXP poly_hmm_est_CPP(SEXP m_R,
                                 SEXP n_mar_R,
                                 SEXP n_ind_R,
                                 SEXP p_R,
                                 SEXP dp_R,
                                 SEXP q_R,
                                 SEXP dq_R,
                                 SEXP g_R,
                                 SEXP rf_R,
                                 SEXP arg_vec_1R,
                                 SEXP loglike_R,
                                 SEXP verbose_R,
                                 SEXP tol_R)
{
  int m = Rcpp::as<int>(m_R);
  int verbose = Rcpp::as<int>(verbose_R);
  int  flag, k1, k2, j, it;
  long double s;
  int n_mar = Rcpp::as<int>(n_mar_R);
  int n_ind = Rcpp::as<int>(n_ind_R);
  long double loglike = Rcpp::as<long double>(loglike_R);
  std::vector<int> p = Rcpp::as<std::vector<int> >(p_R);
  std::vector<int> dp = Rcpp::as<std::vector<int> >(dp_R);
  std::vector<int> q = Rcpp::as<std::vector<int> >(q_R);
  std::vector<int> dq = Rcpp::as<std::vector<int> >(dq_R);
  std::vector<double> g = Rcpp::as<std::vector<double> >(g_R);
  std::vector<double> rf = Rcpp::as<std::vector<double> >(rf_R);
  double tol = Rcpp::as<double>(tol_R);
  std::vector<double> arg_vec_1 = Rcpp::as<std::vector<double> >(arg_vec_1R);
  int gam = nChoosek(m, m/2);
  int gam_pow_2 = gam*gam;
  std::vector<double> rec(gam_pow_2);
  std::vector<double> init(gam_pow_2);
  std::vector<double> a(gam_pow_2);
  std::vector<double> b(gam_pow_2);
  std::vector<double> emit_alpha(gam_pow_2);
  std::vector<double> emit_beta(gam_pow_2);
  std::vector<long double> term(n_ind);
  std::fill(init.begin(), init.end(), 1.0/gam_pow_2);
  std::fill(a.begin(), a.end(), 0);
  std::vector<double> cur_rf(rf.size());
  std::vector<double> pre_calc_prob_alpha((1+m/2)*(1+m/2));
  std::vector<double> pre_calc_prob_beta((1+m/2)*(1+m/2));
  //Dynamic memory allocation using operator new
  long double** alpha = new long double*[gam_pow_2]; 
  for(int i = 0; i < gam_pow_2; ++i){
    alpha[i] = new long double[n_mar];
  }
  long double** beta = new long double*[gam_pow_2];
  for(int i = 0; i < gam_pow_2; ++i){
    beta[i] = new long double[n_mar];
  }
  double maxit=400;

  //caching number of recombinations
  setup_pre_calc_n_rec_cache(m, gam);
  if(verbose==1)
  {
    Rcpp::Rcout << "Ploidy level:" << m << "\n" ;
    Rcpp::Rcout << "Number of individuals:" << n_ind << "\n" ;
  }
  /* begin EM algorithm */
  for(it=0; it<maxit; it++)
  {
    for(j=0; j < n_mar-1; j++)
    {
      cur_rf[j] = rf[j];
      rf[j] = 0.0;
    }
    for(int j=0; j < n_ind; j++)
      term[j]=0;
    for(int ind=0; ind < n_ind; ind++)
    {
      R_CheckUserInterrupt();
      for(int j=0; j < gam_pow_2; j++)
      {
        for(int i=0; i < n_mar; i++)
        {
          alpha[j][i] = 0;
          beta[j][i] = 0;
        }
      }
      //caching emissions for first marker
      emit_alpha = emit_poly(m, ind*(m+1)*n_mar,  dp[0], dp[1], dq[0], dq[1], p, q, g);

      //Initializing alpha and beta
      for(int j=0; j < gam_pow_2; j++)
      {
        alpha[j][0] = init[j] * emit_alpha[j];
        beta[j][n_mar-1] = 1;
      }
      //Forward-backward algorithm
      for(k1=1, k2= n_mar-2; k1 < n_mar; k1++, k2--)
      {
        //caching probabilities
        //necessary info: m, rf
        //for alpha
        for(int i=0; i<(1+m/2); i++)
        {
          for(int j=0; j<(1+m/2); j++)
          {
            pre_calc_prob_alpha[i*(1+m/2)+j]=prob_k1_given_k_lp_lq_m(m,i,j,cur_rf[k1-1]);
          }
        }
        //for beta
        for(int i=0; i<(1+m/2); i++)
        {
          for(int j=0; j<(1+m/2); j++)
          {
            pre_calc_prob_beta[i*(1+m/2)+j]=prob_k1_given_k_lp_lq_m(m,i,j,cur_rf[k2]);
          }
        }
        //caching emissions for marker k+1
        emit_alpha = emit_poly(m, ind*(m+1)*n_mar+k1*(m+1), dp[k1], dp[k1+1], dq[k1], dq[k1+1], p, q, g);
        emit_beta = emit_poly(m, ind*(m+1)*n_mar+(k2+1)*(m+1), dp[k2+1], dp[k2+2], dq[k2+1], dq[k2+2], p, q, g);

        //Induction
        for(int l1 = 0; l1 < gam; l1++)
        {
          for(int l2 = 0; l2 < gam; l2++)
          {
            //alpha
            if(emit_alpha[l1*gam+l2] > TOL)
            {
              a=alpha_ai_dot(pre_calc_prob_alpha, m, gam, gam_pow_2, l1, l2); //
              for(int j=0; j<gam_pow_2; j++)
              {
                alpha[l1*gam+l2][k1] = alpha[l1*gam+l2][k1] + a[j] * alpha[j][k1-1];
              }
              alpha[l1*gam+l2][k1]=alpha[l1*gam+l2][k1]*emit_alpha[l1*gam+l2];
            }
            else
            {
              alpha[l1*gam+l2][k1]=0;
            }
            //beta
            if(emit_beta[l1*gam+l2] > TOL)
            {
              b=alpha_ai_dot(pre_calc_prob_beta, m, gam, gam_pow_2, l1, l2);
              for(int j=0; j<gam_pow_2; j++)
              {
                beta[j][k2] = beta[j][k2] + b[j] * emit_beta[l1*gam+l2] * beta[l1*gam+l2][k2+1];
              }
            }
          }
        }
      }
      for(int k = 0; k < n_mar-1; k++)
      {
        for(int i=0; i<(1+m/2); i++)
        {
          for(int j=0; j<(1+m/2); j++)
          {
            pre_calc_prob_alpha[i*(1+m/2)+j]=prob_k1_given_k_lp_lq_m(m,i,j,cur_rf[k]);
          }
        }
        //caching emissions for marker k+1
        emit_beta = emit_poly(m, ind*(m+1)*n_mar+(k+1)*(m+1), dp[k+1], dp[k+2], dq[k+1], dq[k+2], p, q, g);
        s=0.0;
        for(int l1 = 0; l1 < gam; l1++)
        {
          for(int l2 = 0; l2 < gam; l2++)
          {
            if((emit_beta[l1*gam+l2] > TOL) & (beta[l1*gam+l2][k+1] > TOL))
            {
              a=alpha_ai_dot(pre_calc_prob_alpha, m, gam, gam_pow_2, l1, l2);
              for(int j=0; j<gam_pow_2; j++)
                s = s + alpha[j][k] *  beta[l1*gam+l2][k+1] * emit_beta[l1*gam+l2] * a[j];
            }
          }
        }
        for(int l1 = 0; l1 < gam; l1++)
        {
          for(int l2 = 0; l2 < gam; l2++)
          {
            if((emit_beta[l1*gam+l2] > TOL) & (beta[l1*gam+l2][k+1] > TOL))
            {
              rec=nrec(m, gam, gam_pow_2, l1, l2);
              a=alpha_ai_dot(pre_calc_prob_alpha, m, gam, gam_pow_2, l1, l2);
              for(int j=0; j<gam_pow_2; j++)
                rf[k] +=  (rec[j] * (alpha[j][k] *
                  beta[l1*gam+l2][k+1] * emit_beta[l1*gam+l2] *a[j])/s);
            }
          }
        }
      }
      //Termination
      for(int j=0; j < gam_pow_2; j++)
      {
        term[ind] = term[ind] + alpha[j][n_mar-1];
      }
    }/* loop over individuals */
  if(tol>=0.99) /*this is used to provide the likelihood given a vector of recombination fractions*/
  {
    for(j=0; j < n_mar-1; j++)
    {
      rf[j] = cur_rf[j];
    }
    for(int j=0; j < n_ind; j++)
    {
      loglike += log10(term[j]);
    }
    List z  = List::create(loglike, rf);
    return z ;
  }
  /* rescale */
  for(int j=0; j<n_mar-1; j++)
  {
    rf[j] /= (double)n_ind;
    if(rf[j] < tol/100.0) rf[j] = tol/100.0;
    else if(rf[j] > 0.5-tol/100.0) rf[j] = 0.5-tol/100.0;
  }

  /* check convergence */
  flag=0;
  for(int j=0; j<n_mar-1; j++)
  {
    if(fabs(rf[j] - cur_rf[j]) > tol*(cur_rf[j]+tol*100.0))
    {
      flag = 1;
      break;
    }
  }
  for(int j=0; j<n_mar-1; j++)
  {
    if(cur_rf[j] >= .45)
    {
      flag = 0;
      break;
    }
  }
  if(verbose)
  {
    Rcpp::Rcout << "\t\n";
    for(int j=0; j<n_mar-1; j++)
    {
      Rcpp::Rcout.precision(3);
      Rcpp::Rcout << std::fixed << rf[j] << " ";
    }
  }
  if(!flag) break;
  } /* end EM algorithm */
  /*Loglike computation*/

  for(int j=0; j < n_ind; j++)
  {
    loglike += log10(term[j]);
  }
  // Cleaning up memory
  for(int i = 0; i < gam_pow_2; ++i) {
    delete [] alpha[i];
    delete [] beta[i];
  }
  delete [] alpha;
  List z  = List::create(loglike, rf);
  return z ;
}

RcppExport SEXP calc_genoprob_prior(SEXP m_R,
                                    SEXP n_mar_R,
                                    SEXP n_ind_R,
                                    SEXP p_R,
                                    SEXP dp_R,
                                    SEXP q_R,
                                    SEXP dq_R,
                                    SEXP g_R,
                                    SEXP rf_R,
                                    SEXP probsR,
                                    SEXP loglike_R,
                                    SEXP verbose_R)
{
  int m = Rcpp::as<int>(m_R);
  int verbose = Rcpp::as<int>(verbose_R);
  std::vector<int> p = Rcpp::as<std::vector<int> >(p_R);
  std::vector<int> dp = Rcpp::as<std::vector<int> >(dp_R);
  std::vector<int> q = Rcpp::as<std::vector<int> >(q_R);
  std::vector<int> dq = Rcpp::as<std::vector<int> >(dq_R);
  std::vector<double> g = Rcpp::as<std::vector<double> >(g_R);
  std::vector<double> cur_rf = Rcpp::as<std::vector<double> >(rf_R);
  std::vector<long double> probs = Rcpp::as<std::vector<long double> >(probsR);
  int gam = nChoosek(m, m/2);
  int gam_pow_2 = gam*gam;
  int  k1, k2;
  int n_mar = Rcpp::as<int>(n_mar_R);
  int n_ind = Rcpp::as<int>(n_ind_R);
  std::vector<double> rec(gam_pow_2);
  std::vector<double> init(gam_pow_2);
  std::vector<double> a(gam_pow_2);
  std::vector<double> b(gam_pow_2);
  std::vector<double> emit_alpha(gam_pow_2);
  std::vector<double> emit_beta(gam_pow_2);
  std::fill(init.begin(), init.end(), 1.0/gam_pow_2);
  std::fill(a.begin(), a.end(), 0);
  long double** alpha = new long double*[gam_pow_2]; 
  for(int i = 0; i < gam_pow_2; ++i){
    alpha[i] = new long double[n_mar];
  }
  long double** beta = new long double*[gam_pow_2];
  for(int i = 0; i < gam_pow_2; ++i){
    beta[i] = new long double[n_mar];
  }
  std::vector<double> pre_calc_prob_alpha((1+m/2)*(1+m/2));
  std::vector<double> pre_calc_prob_beta((1+m/2)*(1+m/2));
  int count = 0;

  //caching number of recombinations
  setup_pre_calc_n_rec_cache(m, gam);
  if(verbose==1)
  {
    Rcpp::Rcout << "Ploidy level:" << m << "\n" ;
    Rcpp::Rcout << "Number of individuals:" << n_ind << "\n" ;
  }
  if(verbose) Rcpp::Rcout << "\t";
  for(int ind=0; ind < n_ind; ind++)
  {
    if(verbose){
      Rcpp::Rcout << ".";
      if((ind+1)%50 == 0) Rcpp::Rcout << "\n\t";
    }
    R_CheckUserInterrupt();
    for(int j=0; j < gam_pow_2; j++)
    {
      for(int i=0; i < n_mar; i++)
      {
        alpha[j][i] = 0;
        beta[j][i] = 0;
      }
    }
    //caching emissions for first marker
    emit_alpha = emit_poly(m, ind*(m+1)*n_mar,  dp[0], dp[1], dq[0], dq[1], p, q, g);

    //Initializing alpha and beta
    for(int j=0; j < gam_pow_2; j++)
    {
      alpha[j][0] = init[j] * emit_alpha[j];
      beta[j][n_mar-1] = 1;
    }
    //Forward-backward algorithm
    for(k1=1, k2= n_mar-2; k1 < n_mar; k1++, k2--)
    {
      //caching probabilities
      //necessary info: m, rf
      //for alpha
      for(int i=0; i<(1+m/2); i++)
      {
        for(int j=0; j<(1+m/2); j++)
        {
          pre_calc_prob_alpha[i*(1+m/2)+j]=prob_k1_given_k_lp_lq_m(m,i,j,cur_rf[k1-1]);
        }
      }
      //for beta
      for(int i=0; i<(1+m/2); i++)
      {
        for(int j=0; j<(1+m/2); j++)
        {
          pre_calc_prob_beta[i*(1+m/2)+j]=prob_k1_given_k_lp_lq_m(m,i,j,cur_rf[k2]);
        }
      }
      //caching emissions for marker k+1
      emit_alpha = emit_poly(m, ind*(m+1)*n_mar+k1*(m+1), dp[k1], dp[k1+1], dq[k1], dq[k1+1], p, q, g);
      emit_beta = emit_poly(m, ind*(m+1)*n_mar+(k2+1)*(m+1), dp[k2+1], dp[k2+2], dq[k2+1], dq[k2+2], p, q, g);
      //Induction
      for(int l1 = 0; l1 < gam; l1++)
      {
        for(int l2 = 0; l2 < gam; l2++)
        {
          //alpha
          if(emit_alpha[l1*gam+l2] > TOL)
          {
            a=alpha_ai_dot(pre_calc_prob_alpha, m, gam, gam_pow_2, l1, l2); //
            for(int j=0; j<gam_pow_2; j++)
            {
              alpha[l1*gam+l2][k1] = alpha[l1*gam+l2][k1] + a[j] * alpha[j][k1-1];
            }
            alpha[l1*gam+l2][k1]=alpha[l1*gam+l2][k1]*emit_alpha[l1*gam+l2];
          }
          else
          {
            alpha[l1*gam+l2][k1]=0;
          }
          //beta
          if(emit_beta[l1*gam+l2] > TOL)
          {
            b=alpha_ai_dot(pre_calc_prob_beta, m, gam, gam_pow_2, l1, l2);
            for(int j=0; j<gam_pow_2; j++)
            {
              beta[j][k2] = beta[j][k2] + b[j] * emit_beta[l1*gam+l2] * beta[l1*gam+l2][k2+1];
            }
          }
        }
      }
    }
    for(int k=0; k < n_mar; k++)
    {
      long double w = 0.0;
      for(int j=0; j < gam_pow_2; j++)
        w += alpha[j][k]*beta[j][k];
      for(int j=0; j < gam_pow_2; j++)
      {
        probs[count] = alpha[j][k]*beta[j][k]/w;
        count++;
      }
    }
  }/* loop over individuals */
  // Cleaning up memory
  for(int i = 0; i < gam_pow_2; ++i) {
    delete [] alpha[i];
    delete [] beta[i];
  }
  delete [] alpha;
  List z  = List::create(probs);
  return z ;
}
