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
 File: hmm_elements.cpp 
 
 Description: This file contains the functions used 
 in a hidden Markov model (HMM).
 
 Functions Written by Marcelo Mollinari.
 
 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
 First version: Dec 19, 2013
 Last update: Feb 10, 2018
 */

#include <algorithm>
#include <iostream>
#include <vector>
#include "combinatorial.h"
#include "hmm_elements.h"
#include <math.h>
#include <Rmath.h>
#include <Rcpp.h>
#include <R_ext/PrtUtil.h>

using namespace std;
using namespace Rcpp;

/* FUNCTION: emit_poly
 
 Emission distribution (equation 8 on the paper):
 Computes the probability of observing a molecular phenotype
 given the multialleleic genotype, i.e. P(Observation|Genotype)
 
 m: ploidy
 
 cte: Argument to get the adequate position  on the molecular phenotype (g).
 It aims to position the function in one of the columns of the list "g"
 presented bellow.
 
 ip_k: start position for reading the genotypes of the parent P
 
 ip_k1: end position for reading the genotypes of the parent P
 
 p: Genotype of parent P, i.e. which homologous have the allele in P.
 For example, if an octaploid P has 4 alleles in homologous 1,3,6,7,
 as shown above,
 ---A---
 ---a---
 ---A---
 ---a---
 ---a---
 ---A---
 ---A---
 ---a---
 
 p[ip_k:ip_k1] =(1,3,6,7)
 
 iq_k, iq_k1, q: Analogously for parent Q
 
 g: probability vector for marker genotypes (doses). For instance,
 octaploid, 3 markers, 2 individuals:
 
 $Ind1
 M1  M2  M3
 0 0.0 0.9 0.8
 1 0.8 0.1 0.1
 2 0.2 0.0 0.1
 3 0.0 0.0 0.0
 4 0.0 0.0 0.0
 5 0.0 0.0 0.0
 6 0.0 0.0 0.0
 7 0.0 0.0 0.0
 8 0.0 0.0 0.0
 
 $Ind2
 M1  M2  M3
 0 0.0 0.9 0.8
 1 0.8 0.1 0.1
 2 0.2 0.0 0.1
 3 0.0 0.0 0.0
 4 0.0 0.0 0.0
 5 0.0 0.0 0.0
 6 0.0 0.0 0.0
 7 0.0 0.0 0.0
 8 0.0 0.0 0.0
 
 In this case, g would be:
 [1] 0.0 0.8 0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.9 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.8
 [20] 0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.8 0.2 0.0 0.0 0.0 0.0 0.0 0.0 0.9 0.1
 [39] 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.8 0.1 0.1 0.0 0.0 0.0 0.0 0.0 0.0
 
 ----------------------------------------------------- */

std::vector <double>  emit_poly(int m, int cte, int ip_k, int ip_k1,
                                int iq_k, int iq_k1,
                                std::vector<int>& p,
                                std::vector<int>& q,
                                std::vector<double>& g)
{
  int h = 0, i, j, k, l;
  int gam = nChoosek(m, m/2);
  int gam_2 = gam*gam;
  //double r = 0.0, dkP=ip_k1-ip_k, dkQ=iq_k1-iq_k;
  std::vector<bool> vp(m), vq(m);
  std::vector<double> emit_prob(gam_2);
  std::fill(vp.begin(), vp.end()-m/2, true);
  std::fill(vq.begin(), vq.end()-m/2, true);
  std::vector<double> S(m+1);
  std::fill(S.begin(), S.end(), 0.0);
  std::fill(emit_prob.begin(), emit_prob.end(), 0.0);
  do
  {
    k=0;
    if(p[ip_k]!=0)
    {
      for(i=ip_k; i < ip_k1; i++)
        k=k+vp[p[i]-1];
    }
    do
    {
      l=0;
      if(q[iq_k]!=0)
      {
        for(j=iq_k; j < iq_k1; j++)
          l=l+vq[q[j]-1];
      }
      if(g[l+k+cte]!=0)
        emit_prob[h]=g[l+k+cte];
      h++;
    }
    while (std::prev_permutation(vq.begin(), vq.end()));
  }
  while (std::prev_permutation(vp.begin(), vp.end()));
  return emit_prob;
}

/* FUNCTION: prob_k1_given_k_l_m
 This is equation 5 on the paper
 -----------------------------------------------------
 Calculates the genotypic transition probability based
 on l, whichdenotes the number of recombinant bivalents
 between loci k and k + 1 in one parent.
 */
double prob_k1_given_k_l_m(int m, int l, double rf)
{
  return ((pow((1-rf),(m/2-l))*pow(rf,l))/nChoosek(m/2, l));
}

/* FUNCTION: rec_num
 -----------------------------------------------------
 Returns a matrix containing the number of recombination
 events between loci k and k + 1 in one parent given the
 ploidy level  m.
 */
std::vector<std::vector<double> > rec_num(int m)
{
  int g = nChoosek(m, m/2);
  std::vector<std::vector<double> > R(g);
  for(int i = 0; (unsigned)i < R.size(); ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      R[i].push_back(n_rec_given_genk_and_k1(m, i+1, j+1)/(double)m);
      //R[i].push_back(n_rec_given_genk_and_k1(m, i+1, j+1)/((double)m/2.0));
    }
  }
  return(R);
}

/* FUNCTION: rec_num_denominator
 -----------------------------------------------------
 Returns a matrix containing the number of recombination
 events between loci k and k + 1 in one parent given the
 ploidy level  m.
 */
std::vector<std::vector<int> > rec_num_no_denominator(int m)
{
  int g = nChoosek(m, m/2);
  std::vector<std::vector<int> > R(g);
  for(int i = 0; (unsigned)i < R.size(); ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      R[i].push_back(n_rec_given_genk_and_k1(m, i+1, j+1));
    }
  }
  return(R);
}

/* FUNCTION: transition
 -----------------------------------------------------
 Returns a transition matrix between loci k and k + 1 in
 one parent i.e. Prop(p_{k+1}|p_k), given the ploidy level
 m and the recombination fraction rf.
 */
std::vector<std::vector<double> > transition(int m, double rf)
{
  int g = nChoosek(m, m/2);
  std::vector<std::vector<double> > T(g);
  for(int i = 0; i < g; ++i)
  {
    for(int j = 0; j < g; ++j)
    {
      T[i].push_back(prob_k1_given_k_l_m(m, n_rec_given_genk_and_k1(m, i+1, j+1), rf));
    }
  }
  return(T);
}

/* FUNCTION:  index_func
 * -----------------------------------------------------
 * This function has a similar purpose as the emission function.
 * It return the indices corresponding to 
 * states that should be visited given two vectors indicating
 * which homologous contain the allelic variant.
 * The result is a vector twice the size of states that should be
 * visited. The first half indicates the indices in the transition
 * space in parent P and the second half indicates the indices in
 * the transition space in parent Q
 */

std::vector<std::vector<int> > index_func(int m,
                                          std::vector<int>& p,
                                          std::vector<int>& q)
{
  int s, ip=0, iq=0;
  //int g = nChoosek(m, m/2);
  std::vector<std::vector<int> > v1(p.size()+q.size()+2);
  std::vector<std::vector<int> > v2(p.size()+q.size()+2);
  std::vector<int> vp(m), vq(m);
  std::fill(vp.begin(), vp.end()-m/2, true);
  std::fill(vq.begin(), vq.end()-m/2, true);
  do
  {
    iq=0;
    do
    {
      s=0;
      for(int j=0; (unsigned)j < p.size(); j++)
        if(p[j]>=0) s=s+vp[p[j]];
        for(int j=0; (unsigned)j < q.size(); j++)
          if(q[j]>=0) s=s+vq[q[j]];
          v1[s].push_back(ip);
          v2[s].push_back(iq);
          v1[v1.size()-1].push_back(ip);
          v2[v2.size()-1].push_back(iq);
          iq++;
    }
    while (std::prev_permutation(vq.begin(), vq.end()));
    ip++;
  }
  while (std::prev_permutation(vp.begin(), vp.end()));
  for(int i=0; (unsigned)i < v1.size(); i++)
  {
    v1[i].insert(v1[i].end(), v2[i].begin(), v2[i].end());
  }
  return(v1);
}

/*FUNCTION: init_poly
 ------------------------------------------------------------------
 Description: This function computes the probability of a certain
 genotype dG in a F1 population, given the genotypes (dosage) of its
 parents (dP and dQ) and the ploidy level (m). This is also know as
 polysomic segregation. It returns the log of the probability
 */
double init_poly(int m, int dP, int dQ, int dG)
{
  int j, i;
  double seg = 0.0;
  for(i=0, j=dG; i <= dG; i++, j--)
    seg += R::dhyper(i, dP, m-dP, m/2, 0) * R::dhyper(j, dQ, m-dQ, m/2, 0);
  return(seg);
}

/* FUNCTION: forward
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward(int m,
                            std::vector<double>& fk,
                            std::vector<int>& ik,
                            std::vector<int>& ik1,
                            std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk1);
}
/* FUNCTION: backward
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward(int m,
                             std::vector<double>& fk1,
                             std::vector<int>& ik,
                             std::vector<int>& ik1,
                             std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk);
}



/* FUNCTION: forward_emit (with both informative parents)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit(int m,
                                 std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}
/* FUNCTION: backward (with both informative parents)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit(int m,
                                  std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]] * emit[k1]; 
    }
  }
  return(fk);
}

/* FUNCTION: forward (with one informative parent)
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 */
std::vector<double> forward_emit_one_parent(int m,
                                            std::vector<double>& fk,
                                            std::vector<int>& ik,
                                            std::vector<int>& ik1,
                                            std::vector<double>& emit,
                                            std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}
/* FUNCTION: backward (with one informative parent)
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 */
std::vector<double> backward_emit_one_parent(int m,
                                             std::vector<double>& fk1,
                                             std::vector<int>& ik,
                                             std::vector<int>& ik1,
                                             std::vector<double>& emit,
                                             std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size();
  int ngenk1 = ik1.size();
  std::vector<double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * emit[k1]; 
    }
  }
  return(fk);
}

/* FUNCTION: forward
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 This is using high precision long double variable.
 */
std::vector<long double> forward_highprec(int m,
                                          std::vector<long double>& fk,
                                          std::vector<int>& ik,
                                          std::vector<int>& ik1,
                                          std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk1);
}
/* FUNCTION: backward
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 This is using high precision long double variable.
 */
std::vector<long double> backward_highprec(int m,
                                           std::vector<long double>& fk1,
                                           std::vector<int>& ik,
                                           std::vector<int>& ik1,
                                           std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
  }
  return(fk);
}


/* FUNCTION: forward_emit_highprec
 -----------------------------------------------------
 Classical forward equation presented in Rabiner 1989.
 Multiallelic implementation
 Using high precision long double variable.
 */
std::vector<long double> forward_emit_highprec(int m,
                                               std::vector<long double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk1(ngenk1);
  std::fill(fk1.begin(), fk1.end(), 0.0);
  for(int k1 = 0; k1 < ngenk1; k1++ )
  {
    for(int k = 0; k < ngenk; k++ )
    {
      fk1[k1] = fk1[k1] + fk[k] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]];
    }
    fk1[k1] = fk1[k1] * emit[k1];
  }
  return(fk1);
}
/* FUNCTION: backward_emit_highprec
 -----------------------------------------------------
 Classical backward equation presented in Rabiner 1989.
 Multiallelic implementation
 Using high precision long double variable.
 */
std::vector<long double> backward_emit_highprec(int m,
                                                std::vector<long double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T)
{
  int ngenk = ik.size()/2;
  int ngenk1 = ik1.size()/2;
  std::vector<long double> fk(ngenk);
  std::fill(fk.begin(), fk.end(), 0.0);
  for(int k = 0; k < ngenk; k++ )
  {
    for(int k1 = 0; k1 < ngenk1; k1++ )
    {
      fk[k] =  fk[k] + fk1[k1] * T[ik[k]][ik1[k1]] * T[ik[k+ngenk]][ik1[k1+ngenk1]] * emit[k1]; 
    }
  }
  return(fk);
}
