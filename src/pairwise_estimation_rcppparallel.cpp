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

// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>
#include <float.h> /* DBL_EPSILON */
#include <R_ext/Applic.h>
#include <algorithm>
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include "combinatorial.h"
#include "hmm_elements.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>
#define THRESH 200.0

using namespace std;
using namespace Rcpp;
using namespace RcppParallel;

template<typename T>
vector<T> s(vector<T> const &v, int m, int n) {
   auto first = v.begin() + m;
   auto last = v.begin() + n + 1;
   vector<T> vector(first, last);
   return vector;
}

struct JsDistance : public Worker {
   
  // input matrix to read from
  const RMatrix<double> mrk_pairs;
  const RMatrix<double> geno;
  const std::vector<int> dP;
  const std::vector<int> dQ;
  const std::vector<double> count_vector;
  const std::vector<std::string> count_matrix_rownames;
  const std::vector<int> count_matrix_number;
  const std::vector<int> count_matrix_pos;
  const std::vector<int> count_matrix_length;
  const int m;
  const double tol;
  const int n_ind;
  // output matrix to write to
  RVector<double> out;
   
  // initialize from Rcpp input and output matrixes (the RMatrix class
  // can be automatically converted to from the Rcpp matrix type)
  JsDistance(const NumericMatrix mrk_pairs,
             const NumericMatrix geno,
             const std::vector<int> dP,
             const std::vector<int> dQ,
             const std::vector<double> count_vector,
             const std::vector<std::string> count_matrix_rownames,
             const std::vector<int> count_matrix_number,
             const std::vector<int> count_matrix_pos,
             const std::vector<int> count_matrix_length,
             const int m,
             const double tol,
             const int n_ind,
             NumericVector out)
  : mrk_pairs(mrk_pairs), geno(geno), dP(dP), dQ(dQ), count_vector(count_vector), count_matrix_rownames(count_matrix_rownames), count_matrix_number(count_matrix_number), count_matrix_pos(count_matrix_pos), count_matrix_length(count_matrix_length), m(m), tol(tol), n_ind(n_ind), out(out) {}
   
  // function call operator that work for the specified range (begin/end)
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t k = begin; k < end; k++) {
      int id = (m+1)*(m+1)*(m+1)*dQ[mrk_pairs(1,k)]+(m+1)*(m+1)*dQ[mrk_pairs(0,k)]+(m+1)*dP[mrk_pairs(1,k)]+dP[mrk_pairs(0,k)]+1;
      double min_lod = 1000000.0;
      double min_lod_rf = -1.0;
      double res1=0.0, res2 = 0.0;

      if(count_matrix_number[(id-1)] > 1) //MODIFIED
        {
          RMatrix<double>::Row gen_1 = geno.row(mrk_pairs(0,k));
          RMatrix<double>::Row gen_2 = geno.row(mrk_pairs(1,k));
          // NumericMatrix res(3, count_matrix_number[(id-1)]); //MODIFIED

          // Getting rownames
          std::string mystring = count_matrix_rownames[(id-1)];
          char split = '/';
          std::vector<std::string> z;
          int start = 0;
          int rn = 1;
          // Looping through string
          for(rn = 1; rn < mystring.size(); rn++){
            if( mystring[rn] == split){
              std::string temp = mystring.substr(start, (rn - start));
              z.push_back(temp);
              start = rn+1;
              rn = rn+1;
            }
          }
          std::string temp = mystring.substr(start, (rn - start));
          z.push_back(temp);
          // End of getting rownames

          // Getting length of actual submatrices
          int len_ac = count_matrix_length[(id-1)] / count_matrix_number[(id-1)];

          for(int i=0; i < count_matrix_number[(id-1)]; i++) //MODIFIED
            {
              std::vector<double> count_mat = s(count_vector,
                                                ((count_matrix_pos[(id-1)]-1) + (i*(len_ac))),
                                                ((count_matrix_pos[(id-1)]-1) + ((i+1)*len_ac - 1)));
              std::vector<int> dk(z.size()), dk1(z.size());
              std::string delimiter = " ";
              for(int j=0; j < z.size(); j++)
                {
                  std::string lnames = z[j];
                  dk[j] = std::stoi(lnames.substr(0,lnames.find(delimiter)));
                  dk1[j] = std::stoi(lnames.substr(lnames.find(delimiter)+1, lnames.length()));
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

              int count, count2=0;
              double temp=0.0;
              std::vector<double> Tr((m+2) * (m+2));
              std::fill(Tr.begin(), Tr.end(), 1);
              for(int i = 0; i < dk.size(); i++){
                //Rcpp::Rcout << dk(i) << " " << dk1(i) <<  std::endl;;
                count=0;
                Tr[dk[i] + (dk1[i] * (m+2))]=0.0;
                for(int lp = 0; lp <= m/2; lp++){
                  for(int lq = lp; lq <= m/2; lq++) {
                    Tr[dk[i] + (dk1[i] * (m+2))]=Tr[dk[i] + (dk1[i] * (m+2))] +
                      count_mat[count2 + (count * z.size())] *
                      pow(x,(lq+lp)) *
                      pow((1-x),(m-lq-lp))/
                      (nChoosek(m/2,lq) *
                       nChoosek(m/2,lp));
                    count++;
                  }
                }
                count2++;
              }
              for(int i = 0; i < n_ind; i++)
                temp += log10(Tr[gen_1[i] + (gen_2[i] * (m+2))]);
              fx = -temp ;
              
              // fx = twopt_likelihood_dosage_rcpp(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
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

                count=0;
                count2=0;
                temp=0.0;
                std::fill(Tr.begin(), Tr.end(), 1);
                for(int i = 0; i < dk.size(); i++){
                  //Rcpp::Rcout << dk(i) << " " << dk1(i) <<  std::endl;;
                  count=0;
                  Tr[dk[i] + (dk1[i] * (m+2))]=0;
                  for(int lp = 0; lp <= m/2; lp++){
                    for(int lq = lp; lq <= m/2; lq++) {
                      Tr[dk[i] + (dk1[i] * (m+2))]=Tr[dk[i] + (dk1[i] * (m+2))] +
                        count_mat[count2 + (count * z.size())] *
                        pow(u,(lq+lp)) *
                        pow((1-u),(m-lq-lp))/
                        (nChoosek(m/2,lq) *
                         nChoosek(m/2,lp));
                      count++;
                    }
                  }
                  count2++;
                }  
                for(int i = 0; i < n_ind; i++)
                  temp += log10(Tr[gen_1[i] + (gen_2[i] * (m+2))]);
                fu = -temp ;
              
                // fu = twopt_likelihood_dosage_rcpp(u, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
          
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

              count=0,
              count2=0;
              temp=0.0;
              std::fill(Tr.begin(), Tr.end(), 1);
              for(int i = 0; i < dk.size(); i++){
                //Rcpp::Rcout << dk(i) << " " << dk1(i) <<  std::endl;;
                count=0;
                Tr[dk[i] + (dk1[i] * (m+2))]=0;
                for(int lp = 0; lp <= m/2; lp++){
                  for(int lq = lp; lq <= m/2; lq++) {
                    Tr[dk[i] + (dk1[i] * (m+2))]=Tr[dk[i] + (dk1[i] * (m+2))] +
                      count_mat[count2 + (count * z.size())] *
                      pow(x,(lq+lp)) *
                      pow((1-x),(m-lq-lp))/
                      (nChoosek(m/2,lq) *
                       nChoosek(m/2,lp));
                    count++;
                  }
                }
                count2++;
              }  
              for(int i = 0; i < n_ind; i++)
                temp += log10(Tr[gen_1[i] + (gen_2[i] * (m+2))]);
              res2 = -temp ;
              
              // res2 = twopt_likelihood_dosage_rcpp(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
              if (min_lod > res2){
                min_lod = res2;
                min_lod_rf = x;
              }
            }      
          out[k] = min_lod_rf;
        }
      else
        {
          out[k] = min_lod_rf;
        }
    }       
  }
};


// [[Rcpp::export]]
RcppExport SEXP pairwise_rf_estimation_disc_rcpp(SEXP mrk_pairs_R,
                                                 SEXP m_R,
                                                 SEXP geno_R,
                                                 SEXP dP_R,
                                                 SEXP dQ_R,
                                                 SEXP count_vector_R,
                                                 SEXP count_matrix_rownames_R,
                                                 SEXP count_matrix_number_R,
                                                 SEXP count_matrix_pos_R,
                                                 SEXP count_matrix_length_R,
                                                 SEXP tol_R,
                                                 SEXP threads_R)
{
  NumericMatrix mrk_pairs = mrk_pairs_R;
  NumericMatrix geno = geno_R;
  std::vector<int> dP = as<std::vector<int>>(dP_R);
  std::vector<int> dQ = as<std::vector<int>>(dQ_R);
  std::vector<double> count_vector = as<std::vector<double>>(count_vector_R);
  std::vector<std::string> count_matrix_rownames = as<std::vector<std::string>>(count_matrix_rownames_R);
  std::vector<int> count_matrix_number = as<std::vector<int>>(count_matrix_number_R);
  std::vector<int> count_matrix_pos = as<std::vector<int>>(count_matrix_pos_R);
  std::vector<int> count_matrix_length = as<std::vector<int>>(count_matrix_length_R);
  NumericVector out(mrk_pairs.ncol());
  int m = as<int>(m_R);
  double tol = as<double>(tol_R);
  int threads = as<int>(threads_R);
  int n_ind = geno.ncol();

  JsDistance jsDistance(mrk_pairs, geno, dP, dQ, count_vector, count_matrix_rownames, count_matrix_number, count_matrix_pos, count_matrix_length, m, tol, n_ind, out);
  
  // call it with parallelFor
  parallelFor(0, mrk_pairs.ncol(), jsDistance, 3000);
  
  return out;
}


// // [[Rcpp::export]]
// RcppExport SEXP pairwise_rf_estimation_disc_rcpp(SEXP mrk_pairs_R,
//                                                  SEXP m_R,
//                                                  SEXP geno_R,
//                                                  SEXP dP_R,
//                                                  SEXP dQ_R,
//                                                  SEXP count_vector_R,
//                                                  SEXP count_matrix_rownames_R,
//                                                  SEXP count_matrix_number_R,
//                                                  SEXP count_matrix_pos_R,
//                                                  SEXP count_matrix_length_R,
//                                                  SEXP tol_R)
// {
//   NumericMatrix mrk_pairs = mrk_pairs_R;
//   NumericMatrix geno = geno_R;
//   std::vector<int> dP = as<std::vector<int>>(dP_R);
//   std::vector<int> dQ = as<std::vector<int>>(dQ_R);
//   std::vector<double> count_vector = as<std::vector<double>>(count_vector_R);
//   std::vector<std::string> count_matrix_rownames = as<std::vector<std::string>>(count_matrix_rownames_R);
//   std::vector<int> count_matrix_number = as<std::vector<int>>(count_matrix_number_R);
//   std::vector<int> count_matrix_pos = as<std::vector<int>>(count_matrix_pos_R);
//   std::vector<int> count_matrix_length = as<std::vector<int>>(count_matrix_length_R);
//   NumericMatrix out(geno.nrow(), geno.nrow());
//   int m = as<int>(m_R);
//   double tol = as<double>(tol_R);
//   int n_ind = geno.ncol();
//   for(int k=0; k < mrk_pairs.ncol(); k++)
//   {
//     int id = (m+1)*(m+1)*(m+1)*dQ[mrk_pairs(1,k)]+(m+1)*(m+1)*dQ[mrk_pairs(0,k)]+(m+1)*dP[mrk_pairs(1,k)]+dP[mrk_pairs(0,k)]+1;
//     double min_lod = 1000000.0;
//     double min_lod_rf = -1.0;

//     if(count_matrix_number[(id-1)] > 1) //MODIFIED
//     {
//       NumericVector gen_1 = geno( mrk_pairs(0,k), _);
//       NumericVector gen_2 = geno( mrk_pairs(1,k), _);
//       NumericMatrix res(3, count_matrix_number[(id-1)]); //MODIFIED
//       // NumericMatrix res(3, temp_list.size()); //MODIFIED
//       // CharacterVector zn = temp_list.attr( "names" ) ; //MODIFIED
//       // colnames(res)=zn; //MODIFIED

//       // Getting rownames
//       std::string mystring = count_matrix_rownames[(id-1)];
//       char split = '/';
//       std::vector<std::string> z;
//       int start = 0;
//       int rn = 1;
//       // Looping through string
//       for(rn = 1; rn < mystring.size(); rn++){
//         if( mystring[rn] == split){
//           std::string temp = mystring.substr(start, (rn - start));
//           z.push_back(temp);
//           start = rn+1;
//           rn = rn+1;
//         }
//       }
//       std::string temp = mystring.substr(start, (rn - start));
//       z.push_back(temp);
//       // End of getting rownames

//       // Getting length of actual submatrices
//       int len_ac = count_matrix_length[(id-1)] / count_matrix_number[(id-1)];

//       for(int i=0; i < count_matrix_number[(id-1)]; i++) //MODIFIED
//       {
//         std::vector<double> count_mat = s(count_vector,
//                                           ((count_matrix_pos[(id-1)]-1) + (i*(len_ac))),
//                                           ((count_matrix_pos[(id-1)]-1) + ((i+1)*len_ac - 1)));
//         std::vector<int> dk(z.size()), dk1(z.size());
//         std::string delimiter = " ";
//         for(int j=0; j < z.size(); j++)
//         {
//           std::string lnames = z[j];
//           dk[j] = std::stoi(lnames.substr(0,lnames.find(delimiter)));
//           dk1[j] = std::stoi(lnames.substr(lnames.find(delimiter)+1, lnames.length()));
//         }
//         //an approximation  x  to the point where  f  attains a minimum  on
//         //the interval  (a,b)  is determined.
//         // Adapted from Brent_fmin function, which can be found in R/src/library/stats/src/optimize.c
//         //   Copyright (C) 1995, 1996  Robert Gentleman and Ross Ihaka
//         //   Copyright (C) 2003-2004  The R Foundation
//         //   Copyright (C) 1998--2014-2018 The R Core Team
//         // This function subprogram is a slightly modified  version  of  the
//         // Algol  60 procedure  localmin  given in Richard Brent, Algorithms for
//         // Minimization without Derivatives, Prentice-Hall, Inc. (1973).
//         // Brent's Minimization Procedure starts
//         //  c is the squared inverse of the golden ratio
//         const double c = (3. - sqrt(5.)) * .5;
//         // Local variables
//         double a, b, d, e, p, q, r, u, v, w, x;
//         double t2, fu, fv, fw, fx, xm, eps, tol1, tol3;
//         //  eps is approximately the square root of the relative machine precision.
//         eps = DBL_EPSILON;
//         tol1 = eps + 1.;// the smallest 1.000... > 1
//         eps = sqrt(eps);
        
//         a = 0.0;
//         b = 0.5;
//         v = a + c * (b - a);
//         w = v;
//         x = v;
        
//         d = 0.;// -Wall
//         e = 0.;
//         fx = twopt_likelihood_dosage_rcpp(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
//         fv = fx;
//         fw = fx;
//         tol3 = tol / 3.;
        
//         //  main loop starts here -----------------------------------
        
//         for(;;) {
//           xm = (a + b) * .5;
//           tol1 = eps * fabs(x) + tol3;
//           t2 = tol1 * 2.;
          
//           // check stopping criterion
          
//           if (fabs(x - xm) <= t2 - (b - a) * .5) break;
//           p = 0.;
//           q = 0.;
//           r = 0.;
//           if (fabs(e) > tol1) { // fit parabola
            
//             r = (x - w) * (fx - fv);
//             q = (x - v) * (fx - fw);
//             p = (x - v) * q - (x - w) * r;
//             q = (q - r) * 2.;
//             if (q > 0.) p = -p; else q = -q;
//             r = e;
//             e = d;
//           }
          
//           if (fabs(p) >= fabs(q * .5 * r) ||
//               p <= q * (a - x) || p >= q * (b - x)) { // a golden-section step
            
//             if (x < xm) e = b - x; else e = a - x;
//             d = c * e;
//           }
//           else { // a parabolic-interpolation step
            
//             d = p / q;
//             u = x + d;
            
//             // f must not be evaluated too close to ax or bx
            
//             if (u - a < t2 || b - u < t2) {
//               d = tol1;
//               if (x >= xm) d = -d;
//             }
//           }
          
//           // f must not be evaluated too close to x
          
//           if (fabs(d) >= tol1)
//             u = x + d;
//           else if (d > 0.)
//             u = x + tol1;
//           else
//             u = x - tol1;
          
//           fu = twopt_likelihood_dosage_rcpp(u, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
          
//           //  update  a, b, v, w, and x
          
//           if (fu <= fx) {
//             if (u < x) b = x; else a = x;
//             v = w;    w = x;   x = u;
//             fv = fw; fw = fx; fx = fu;
//           } else {
//             if (u < x) a = u; else b = u;
//             if (fu <= fw || w == x) {
//               v = w; fv = fw;
//               w = u; fw = fu;
//             } else if (fu <= fv || v == x || v == w) {
//               v = u; fv = fu;
//             }
//           }
//         }
//         // Brent's Minimization Procedure ends
//         res(0,i) = x;
//         res(1,i) = twopt_likelihood_dosage_rcpp(x, m, n_ind, dP[mrk_pairs(0,k)], dQ[mrk_pairs(0,k)], dk, dk1, gen_1, gen_2, count_mat, z.size());
//         if (min_lod > res(1,i)){
//           min_lod = res(1,i);
//           min_lod_rf = x;
//         }
//       }
//       // out(k)=res;

//       // // Checking for min LOD and getting rf for matrix (fix it)
//       // int pos_max = 0;
//       // for(int pos=0; pos < res.ncol(); pos++){
//       //   double max_lod = 10000.0;
//       //   if(max_lod > res(1,pos)){
//       //     max_lod = res(1,pos);
//       //     pos_max = pos;
//       //   }
//       // }
      
//       out(mrk_pairs(0,k),mrk_pairs(1,k)) = out(mrk_pairs(1,k),mrk_pairs(0,k)) = min_lod_rf;
//     }
//     else
//     {
//       // NumericVector d_out(4);
//       // d_out(0)=dP[mrk_pairs(0,k)];
//       // d_out(1)=dP[mrk_pairs(1,k)];
//       // d_out(2)=dQ[mrk_pairs(0,k)];
//       // d_out(3)=dQ[mrk_pairs(1,k)];
//       // out(k)=d_out;
//       out(mrk_pairs(0,k),mrk_pairs(1,k)) = out(mrk_pairs(1,k),mrk_pairs(0,k)) = min_lod_rf;
//     }
//   }
//   return(out);
// }

