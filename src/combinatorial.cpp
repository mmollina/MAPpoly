/*
  Polymap: a package to construct genetic maps in autopolyploids
  Copyright (C) 2014 Marcelo Mollinari

    This file is part of Polymap.

    Polymap is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
  File: combinatorial_functions.cpp
  Description: Set of combinatorial functions to be used with mappoly

  Functions mostly Written by Marcelo Mollinari.

  Functions allocate_double and allocate_alpha were
  written by Karl W Broman and can be found in the R package qtl
  Copyright (c) 2001-10, Karl W Broman

  Escola Superior de Agricultura "Luiz de Queiroz"
  Departamento de Genética - São Paulo, Brazil
  Contact: mmollina@usp.br
  First version: Dec 19, 2013
  Last update: Jul 31, 2014
*/

#include <R.h>
#include <algorithm>
#include <iostream>
#include <vector>
#include "combinatorial.h"
#include <math.h>
#include <Rmath.h>
#include <R_ext/PrtUtil.h>
#define THRESH 200.0

/**********************************************************************
 *
 * allocate_double
 *
 * Allocate space for a vector of doubles
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_double(int n, double **vector)
{
  *vector = (double *)R_alloc(n, sizeof(double));
}

/**********************************************************************
 *
 * allocate_alpha
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_alpha(int n_pos, int n_gen, double ***alpha)
{
    int i;

    *alpha = (double **)R_alloc(n_gen, sizeof(double *));

    (*alpha)[0] = (double *)R_alloc(n_gen*n_pos, sizeof(double));

    for(i=1; i< n_gen; i++)
        (*alpha)[i] = (*alpha)[i-1] + n_pos;
}

/**********************************************************************
 *
 * allocate_genoprob
 *
 * Allocate space for alpha and beta matrices
 *
 * Afterwards, indexed like Genoprob[gen][mar][ind]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_genoprob_long(int n_pos, int n_gen, int n_ind, long double ****Genoprob)
{
  int i, j;
  long double **a;

  *Genoprob = (long double ***)R_alloc(n_gen, sizeof(long double **));

  a = (long double **)R_alloc(n_pos*n_gen, sizeof(long double *));

  (*Genoprob)[0] = a;
  for(i=1; i< n_gen; i++)
    (*Genoprob)[i] = (*Genoprob)[i-1]+n_pos;
}


/**********************************************************************
 *
 * allocate_alpha_long
 *
 * Allocate space for alpha and beta matrices (long double version)
 *
 * Afterwards, indexed like alpha[gen][mar]
 *
 * Allocation done by R_alloc, so that R does the cleanup.
 *
 **********************************************************************/
void allocate_alpha_long(int n_pos, int n_gen, long double ***alpha)
{
  int i;

  *alpha = (long double **)R_alloc(n_gen, sizeof(long double *));

  (*alpha)[0] = (long double *)R_alloc(n_gen*n_pos, sizeof(long double));

  for(i=1; i< n_gen; i++)
    (*alpha)[i] = (*alpha)[i-1] + n_pos;
}

/* FUNCTION: nChoosek
   -----------------------------------------------------
   The famous binomial coefficient
 */

int nChoosek(int n, int k)
{
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    int result = n;
    for( int i = 2; i <= k ; ++i )
    {
        result *= (n-i+1);
        result /= i;
    }
    return result;
}


/*
  FUNCTION: n_rec_given_genk_and_k1
  -----------------------------------------------------
  Given two boolean vectors representing two genotypes k and k+1, this
  function returns the number of recombinants in a gamete for a
  specific linkage phase. For example, vector [1 1 1 0 0 0] represents
  the genotype P_k^1 P_k^2 and P_k^3 and vector [1 1 0 1 0 0]
  represents the genotype P_{k+1}^1, P_{k+1}^2, P_{k+1}^4.  The number
  of recombinant events in this example is 1.
 */

int n_rec_given_genk_and_k1(int ploidy, int index1, int index2)
{
    int i, result = 0;
    std::vector<bool> vec1(ploidy), vec2(ploidy);
    std::fill(vec1.begin(), vec1.end()-ploidy/2, false);
    std::fill(vec2.begin(), vec2.end()-ploidy/2, false);
    vec1=get_boolean_vec_from_lexicographical_index(ploidy, index1);
    vec2=get_boolean_vec_from_lexicographical_index(ploidy, index2);
    for(i=0; i < ploidy; i++)
    {
        if((vec1[i]+vec2[i]) == 2)
        {
            result++;
        }
    }
    result = ploidy/2 - result;
    return result;
}



/* FUNCTION: prob_k1_given_k_lp_lq_m
   This is equation 6 on the paper
   -----------------------------------------------------
   Calculates the genotypic transition probability based on l_P and
   l_Q l_P and l_Q denote the number of recombinant bivalents between
   loci k and k + 1 in parents P and Q respectively.
 */
double prob_k1_given_k_lp_lq_m(int m,
                               int lp,
                               int lq,
                               double rf)
{
    return ((pow((1-rf),(m/2-lp))*pow(rf,lp))/nChoosek(m/2, lp) *
            (pow((1-rf),(m/2-lq))*pow(rf,lq))/nChoosek(m/2, lq));
}


/* FUNCTION: boolean_lexicographic_k_choose_m_and_collapse This is the
   algorithm 2 on the paper.
   -----------------------------------------------------
   This function combines the adequated conditional probabilities in
   order to make the reduction of diomensionality for the two-point
   analyses. This function returns f(m, lP, lQ, wkP, wQk) presented on
   equation 18. IMPORTANT: Notice that here, the last element on the
   vector 'counts' in a normalization constant.
*/
std::vector <int>  boolean_lexicographic_k_choose_m_and_collapse(int ploidy,
        std::vector<int>& which_homologous_mk1,
        std::vector<int>& which_homologous_mk2,
        int gen_prog_mk1,
        int gen_prog_mk2)
{
    int m = nChoosek(ploidy, ploidy/2);
    int i1, i2, j1 = 0, nrow = 0;
    std::vector<bool> vec1(ploidy), vec2(ploidy);
    std::vector<int> pos1(m), pos2(m);
    std::vector<int> counts(1+ploidy/2);
    std::fill(vec1.begin(), vec1.end()-ploidy/2, true);
    std::fill(pos1.begin(), pos1.end(), 0);
    std::fill(counts.begin(), counts.end(), 0);
    do
    {
        for(i1=0; i1<(int)which_homologous_mk1.size(); i1++)
        {
            pos1[j1] += (int)vec1[which_homologous_mk1[i1]];
        }
        if(gen_prog_mk1 == pos1[j1])
        {
            nrow++;
            std::fill(pos2.begin(), pos2.end(), 0);
            std::fill(vec2.begin(), vec2.end()-ploidy/2, true);
            int j2 = 0;
            do
            {
                for(i2=0; i2<(int)which_homologous_mk2.size(); i2++)
                {
                    pos2[j2] += (int)vec2[which_homologous_mk2[i2]];
                }
                if(gen_prog_mk2 == pos2[j2])
                    counts[n_rec_given_genk_and_k1(ploidy,j1+1,j2+1)]++; //compare strings: much faster
                j2++;
            }
            while (std::prev_permutation(vec2.begin(), vec2.end()));
        }
        j1++;
    }
    while (std::prev_permutation(vec1.begin(), vec1.end()));
    //for(i1=0; i1 < 1+ploidy/2; i1++)
    //  counts[i1] /= nrow;
    counts.push_back(nrow);
    return counts;
}


/* FUNCTION: get_boolean_vec_from_lexicographical_index
   This is algotithm 1 on the paper
   -----------------------------------------------------
   This function takes as arguments the ploidy level and a
   lexicographical index and returns the boolean lexicographical
   combination for that index (in a boolean vector). It is importante
   to notice that the algorithm does not calculate all possible
   lexicographical combinations to get the requested combination.
 */
std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index)
{
    int i, j, increment, sentinel;
    std::vector<bool> vec(ploidy+1);
    i=0;
    j=1;
    increment=0;
    sentinel=0;
    std::fill(vec.begin(), vec.end(), 0);
    while(sentinel < ploidy/2)
    {
        if(index > nChoosek((ploidy-j), (ploidy/2 - (i+1))) + increment)
        {
            vec[j-1]=0;
            increment += nChoosek((ploidy-j), (ploidy/2 - (i+1)));
        }
        else
        {
            vec[j-1]=1;
            i++;
        }
        sentinel += vec[j-1];
        j++;
    }
    return vec;
}
//end of file
