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
 File: two_pts_est.cpp
 
 Functions Written partially by Marcelo Mollinari.
 
 Part of this function was adapted from Brent_fmin function, 
 which can be found in R/src/library/stats/src/optimize.c
 Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
 Copyright (C) 2003-2004  The R Foundation
 Copyright (C) 1998--2014-2018 The R Core Team
 
 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version: Dec 19, 2013
 Last update: Jul 31, 2014
 */

#include <float.h> /* DBL_EPSILON */
#include <R_ext/Applic.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include "combinatorial.h"
#include "hmm_elements.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#define THRESH 200.0

using namespace std;
using namespace Rcpp;

double twopt_likelihood_dosage(double rf, int m, int n_ind,
                               int dP, int dQ,
                               Rcpp::NumericVector dk,
                               Rcpp::NumericVector dk1,
                               Rcpp::NumericVector gen_1,
                               Rcpp::NumericVector gen_2,
                               Rcpp::NumericMatrix count_mat)
{
  int count, count2=0;
  double temp=0.0;
  Rcpp::NumericMatrix Tr(m+2, m+2);
  std::fill(Tr.begin(), Tr.end(), 1);
  for(int i = 0; i < dk.size(); i++){
    //Rcpp::Rcout << dk(i) << " " << dk1(i) <<  std::endl;;
    count=0;
    Tr(dk(i),dk1(i))=0;
    for(int lp = 0; lp <= m/2; lp++){
      for(int lq = lp; lq <= m/2; lq++) {
        Tr(dk(i),dk1(i))=Tr(dk(i),dk1(i)) +
          count_mat(count2, count) *
          pow(rf,(lq+lp)) *
          pow((1-rf),(m-lq-lp))/
            (nChoosek(m/2,lq) *
              nChoosek(m/2,lp));
        count++;
      }
    }
    count2++;
  }
  
  for(int i = 0; i < n_ind; i++)
    temp += log10(Tr(gen_1(i),gen_2(i)));
  return -temp ;
}

double twopt_likelihood_probability(double rf, int m, int n_ind,
                                    int dP, int dQ,
                                    Rcpp::NumericVector dk,
                                    Rcpp::NumericVector dk1,
                                    double** prob_mat1,
                                    double** prob_mat2,
                                    Rcpp::NumericMatrix count_mat)
{
  int count, count2=0;
  double temp=0.0;
  Rcpp::NumericMatrix Tr(m+1, m+1);
  std::fill(Tr.begin(), Tr.end(), 1);
  for(int i = 0; i < dk.size(); i++){
    count=0;
    Tr(dk(i),dk1(i))=0;
    for(int lp = 0; lp <= m/2; lp++){
      for(int lq = lp; lq <= m/2; lq++) {
        Tr(dk(i),dk1(i))=Tr(dk(i),dk1(i)) +
          count_mat(count2, count) *
          pow(rf,(lq+lp)) *
          pow((1-rf),(m-lq-lp))/
            (nChoosek(m/2,lq) *
              nChoosek(m/2,lp));
        count++;
      }
    }
    count2++;
  }
  for(int i = 0; i < n_ind; i++){
    double temp1 = 0;
    for(int j = 0; j <= m; j++){
      double temp2 = 0;
      for(int l = 0; l <= m; l++){
        temp2 += prob_mat1[l][i] * Tr(l,j);
      }
      temp1 += temp2 * prob_mat2[j][i];
    }
    temp += log10(temp1);
  }
  return -temp ;
}

RcppExport SEXP pairwise_rf_estimation_disc(SEXP m_R,
                                            SEXP mrk_pairs_R,
                                            SEXP geno_R,
                                            SEXP dP_R,
                                            SEXP dQ_R,
                                            SEXP count_cache_R,
                                            SEXP tol_R)
{
  Rcpp::NumericMatrix mrk_pairs = Rcpp::as<Rcpp::NumericMatrix>(mrk_pairs_R);
  Rcpp::NumericMatrix geno = Rcpp::as<Rcpp::NumericMatrix>(geno_R);
  Rcpp::NumericVector dP = Rcpp::as<Rcpp::NumericVector>(dP_R);
  Rcpp::NumericVector dQ = Rcpp::as<Rcpp::NumericVector>(dQ_R);
  Rcpp::List count_cache = Rcpp::as<Rcpp::List>(count_cache_R);
  Rcpp::NumericVector d_pair(4);
  Rcpp::List out(mrk_pairs.ncol());
  int m = Rcpp::as<int>(m_R);
  double tol = Rcpp::as<double>(tol_R);
  int n_ind = geno.ncol();
  for(int k=0; k < mrk_pairs.ncol(); k++)
  {
    //Rcpp::Rcout << mrk_pairs(0,k) << " - " << mrk_pairs(1,k) <<  std::endl;
    int id = (m+1)*(m+1)*(m+1)*dQ[mrk_pairs(1,k)]+(m+1)*(m+1)*dQ[mrk_pairs(0,k)]+(m+1)*dP[mrk_pairs(1,k)]+dP[mrk_pairs(0,k)]+1;
    //Rcpp::Rcout << "id: " << id <<  std::endl;
    Rcpp::List temp_list = count_cache[(id-1)];
    //Rcpp::List temp_list = count_cache[k];
    if(temp_list.size() > 1)
    {
      NumericVector gen_1 = geno( mrk_pairs(0,k), _);
      NumericVector gen_2 = geno( mrk_pairs(1,k), _);
      Rcpp::NumericMatrix res(3, temp_list.size());
      Rcpp::CharacterVector zn = temp_list.attr( "names" ) ;
      colnames(res)=zn;
      for(int i=0; i < temp_list.size(); i++)
      {
        Rcpp::NumericMatrix count_mat = temp_list[i] ;
        Rcpp::List dimnames = count_mat.attr( "dimnames" ) ;
        Rcpp::CharacterVector z = dimnames[0];
        Rcpp::NumericVector dk(z.size()), dk1(z.size());
        std::string delimiter = " ";
        for(int j=0; j < z.size(); j++)
        {
          std::string lnames = Rcpp::as<std::string>(z(j));
          dk(j) = std::stoi(lnames.substr(0,lnames.find(delimiter)));
          dk1(j) = std::stoi(lnames.substr(lnames.find(delimiter)+1, lnames.length()));
        }
        //an approximation  x  to the point where  f  attains a minimum  on
        //the interval  (a,b)  is determined.
        // Adapted from Brent_fmin function, which can be found in R/src/library/stats/src/optimize.c
        //   Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
        //   Copyright (C) 2003-2004  The R Foundation
        //   Copyright (C) 1998--2014-2018 The R Core Team
        // This function subprogram is a slightly modified  version  of  the
        // Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
        // Minimization without Derivatives, Prentice-Hall, Inc. (1973).
        // Brent's Minimization Procedure starts
        //  c is the squared inverse of the golden ratio
        const double c = (3. - sqrt(5.)) * .5;
        // Local variables
        double a, b, d, e, p, q, r, u, v, w, x;
        double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
        //  eps is approximately the square root of the relative machine precision.
        eps = DBL_EPSILON;
        tol1 = eps + 1.;// the smallest 1.000... > 1
        eps = sqrt(eps);
        
        a = 0.0;
        b = 0.5;
        v = a + c * (b - a);
        w = v;
        x = v;
        
        d = 0.;// -Wall
        e = 0.;
        fx = twopt_likelihood_dosage(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat);
        //twopt_likelihood (x, m, n_ind, dk, dk1, q_mendel, gen_1, gen_2, count_mat);//(*f)(x, info);
        fv = fx;
        fw = fx;
        tol3 = tol / 3.;
        
        //  main loop starts here -----------------------------------
        
        for(;;) {
          xm = (a + b) * .5;
          tol1 = eps * fabs(x) + tol3;
          t2 = tol1 * 2.;
          
          // check stopping criterion
          
          if (fabs(x - xm) <= t2 - (b - a) * .5) break;
          p = 0.;
          q = 0.;
          r = 0.;
          if (fabs(e) > tol1) { // fit parabola
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
          }
          
          if (fabs(p) >= fabs(q * .5 * r) ||
              p <= q * (a - x) || p >= q * (b - x)) { // a golden-section step
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
          }
          else { // a parabolic-interpolation step
            
            d = p / q;
            u = x + d;
            
            // f must not be evaluated too close to ax or bx
            
            if (u - a < t2 || b - u < t2) {
              d = tol1;
              if (x >= xm) d = -d;
            }
          }
          
          // f must not be evaluated too close to x
          
          if (fabs(d) >= tol1)
            u = x + d;
          else if (d > 0.)
            u = x + tol1;
          else
            u = x - tol1;
          
          fu = twopt_likelihood_dosage(u, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat);
          //twopt_likelihood(u, m, n_ind, dk, dk1, q_mendel, gen_1, gen_2, count_mat);//(*f)(u, info);
          
          //  update  a, b, v, w, and x
          
          if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
          } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
              v = w; fv = fw;
              w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v = u; fv = fu;
            }
          }
        }
        // Brent's Minimization Procedure ends
        res(0,i) = x;
        res(1,i) = twopt_likelihood_dosage(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat);
        res(2,i) = twopt_likelihood_dosage(0.5, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat);
      }
      out(k)=res;
    }
    else
    {
      Rcpp::NumericVector d_out(4);
      d_out(0)=dP[mrk_pairs(0,k)];
      d_out(1)=dP[mrk_pairs(1,k)];
      d_out(2)=dQ[mrk_pairs(0,k)];
      d_out(3)=dQ[mrk_pairs(1,k)];
      out(k)=d_out;
    }
  }
  return(out);
}

RcppExport SEXP pairwise_rf_estimation_prob(SEXP m_R,
                                            SEXP mrk_pairs_R,
                                            SEXP geno_dim_R,
                                            SEXP geno_prob_R,
                                            SEXP dP_R,
                                            SEXP dQ_R,
                                            SEXP count_cache_R,
                                            SEXP tol_R)
{
  Rcpp::NumericMatrix mrk_pairs = Rcpp::as<Rcpp::NumericMatrix>(mrk_pairs_R);
  Rcpp::NumericVector geno_dim = Rcpp::as<Rcpp::NumericVector>(geno_dim_R);
  Rcpp::NumericVector geno_prob = Rcpp::as<Rcpp::NumericVector>(geno_prob_R);
  int n_mrk = geno_dim[0]; //n.mrk
  int n_dos = geno_dim[1]; //m +1
  int n_ind = geno_dim[2]; //n.ind
  int l = 0;
  double ***G = new double**[n_mrk];
  for(int i = 0; i < n_mrk; i++){
    G[i] = new double*[n_dos];
    for(int j = 0; j < n_dos; j++){
      G[i][j] = new double[n_ind];
      for(int k = 0; k < n_ind; k++){
        G[i][j][k] = 0;
      }
    }
  }
  
  for(int k = 0; k < n_ind; k++){
    for(int j = 0; j < n_dos; j++){
      for(int i = 0; i < n_mrk; i++){
        G[i][j][k] = geno_prob[l];
        l++;
      }
    }
  }
  
  Rcpp::NumericVector dP = Rcpp::as<Rcpp::NumericVector>(dP_R);
  Rcpp::NumericVector dQ = Rcpp::as<Rcpp::NumericVector>(dQ_R);
  Rcpp::List count_cache = Rcpp::as<Rcpp::List>(count_cache_R);
  Rcpp::NumericVector d_pair(4);
  Rcpp::List out(mrk_pairs.ncol());
  int m = Rcpp::as<int>(m_R);
  double tol = Rcpp::as<double>(tol_R);
  
  for(int k=0; k < mrk_pairs.ncol(); k++)
  {
    //Rcpp::Rcout << mrk_pairs(0,k) << " - " << mrk_pairs(1,k) <<  std::endl;
    int id = (m+1)*(m+1)*(m+1)*dQ[mrk_pairs(1,k)]+(m+1)*(m+1)*dQ[mrk_pairs(0,k)]+(m+1)*dP[mrk_pairs(1,k)]+dP[mrk_pairs(0,k)]+1;
    //Rcpp::Rcout << "id: " << id <<  std::endl;
    Rcpp::List temp_list = count_cache[(id-1)];
    //Rcpp::List temp_list = count_cache[k];
    if(temp_list.size() > 1)
    {
      int i1 = mrk_pairs(0,k);
      int i2 = mrk_pairs(1,k);
      
      Rcpp::NumericMatrix res(3, temp_list.size());
      Rcpp::CharacterVector zn = temp_list.attr( "names" ) ;
      colnames(res)=zn;
      for(int i=0; i < temp_list.size(); i++)
      {
        Rcpp::NumericMatrix count_mat = temp_list[i] ;
        Rcpp::List dimnames = count_mat.attr( "dimnames" ) ;
        Rcpp::CharacterVector z = dimnames[0];
        Rcpp::NumericVector dk(z.size()), dk1(z.size());
        std::string delimiter = " ";
        for(int j=0; j < z.size(); j++)
        {
          std::string lnames = Rcpp::as<std::string>(z(j));
          dk(j) = std::stoi(lnames.substr(0,lnames.find(delimiter)));
          dk1(j) = std::stoi(lnames.substr(lnames.find(delimiter)+1, lnames.length()));
        }
        //an approximation  x  to the point where  f  attains a minimum  on
        //the interval  (a,b)  is determined.
        // Adapted from Brent_fmin function, which can be found in R/src/library/stats/src/optimize.c
        //   Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
        //   Copyright (C) 2003-2004  The R Foundation
        //   Copyright (C) 1998--2014-2018 The R Core Team
        // This function subprogram is a slightly modified  version  of  the
        // Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
        // Minimization without Derivatives, Prentice-Hall, Inc. (1973).
        // Brent's Minimization Procedure starts
        //  c is the squared inverse of the golden ratio
        const double c = (3. - sqrt(5.)) * .5;
        // Local variables
        double a, b, d, e, p, q, r, u, v, w, x;
        double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
        //  eps is approximately the square root of the relative machine precision.
        eps = DBL_EPSILON;
        tol1 = eps + 1.;// the smallest 1.000... > 1
        eps = sqrt(eps);
        
        a = 0.0;
        b = 0.5;
        v = a + c * (b - a);
        w = v;
        x = v;
        
        d = 0.;// -Wall
        e = 0.;
        
        fx = twopt_likelihood_probability(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, *(G + i1), *(G + i2), count_mat);
        
        fv = fx;
        fw = fx;
        tol3 = tol / 3.;
        
        //  main loop starts here -----------------------------------
        
        for(;;) {
          xm = (a + b) * .5;
          tol1 = eps * fabs(x) + tol3;
          t2 = tol1 * 2.;
          
          // check stopping criterion
          
          if (fabs(x - xm) <= t2 - (b - a) * .5) break;
          p = 0.;
          q = 0.;
          r = 0.;
          if (fabs(e) > tol1) { // fit parabola
            
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = (q - r) * 2.;
            if (q > 0.) p = -p; else q = -q;
            r = e;
            e = d;
          }
          
          if (fabs(p) >= fabs(q * .5 * r) ||
              p <= q * (a - x) || p >= q * (b - x)) { // a golden-section step
            
            if (x < xm) e = b - x; else e = a - x;
            d = c * e;
          }
          else { // a parabolic-interpolation step
            
            d = p / q;
            u = x + d;
            
            // f must not be evaluated too close to ax or bx
            
            if (u - a < t2 || b - u < t2) {
              d = tol1;
              if (x >= xm) d = -d;
            }
          }
          
          // f must not be evaluated too close to x
          
          if (fabs(d) >= tol1)
            u = x + d;
          else if (d > 0.)
            u = x + tol1;
          else
            u = x - tol1;
          
          fu = twopt_likelihood_probability(u, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, *(G + i1), *(G + i2), count_mat);
          
          //  update  a, b, v, w, and x
          
          if (fu <= fx) {
            if (u < x) b = x; else a = x;
            v = w;    w = x;   x = u;
            fv = fw; fw = fx; fx = fu;
          } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
              v = w; fv = fw;
              w = u; fw = fu;
            } else if (fu <= fv || v == x || v == w) {
              v = u; fv = fu;
            }
          }
        }
        // Brent's Minimization Procedure ends
        res(0,i) = x;
        res(1,i) = twopt_likelihood_probability(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, *(G + i1), *(G + i2), count_mat);
        res(2,i) = twopt_likelihood_probability(0.5, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, *(G + i1), *(G + i2), count_mat);
      }
      out(k)=res;
    }
    else
    {
      Rcpp::NumericVector d_out(4);
      d_out(0)=dP[mrk_pairs(0,k)];
      d_out(1)=dP[mrk_pairs(1,k)];
      d_out(2)=dQ[mrk_pairs(0,k)];
      d_out(3)=dQ[mrk_pairs(1,k)];
      out(k)=d_out;
    }
  }
  // Cleaning up memory
  delete [] G;
  return(out);
}
//end of file two_pts_est.cpp
