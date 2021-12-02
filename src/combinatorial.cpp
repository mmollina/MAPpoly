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
  File: combinatorial_functions.cpp
  Description: Set of combinatorial functions to be used with mappoly

  Functions Written by Marcelo Mollinari.

 Bioinformatics Research Center
 Department of Horticultural Science
 North Carolina State University
 Contact: mmollin@ncsu.edu
  First version: Dec 19, 2013
  Last update:   Oct 05, 2021
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
   l_Q, where l_P and l_Q denote the number of recombinant bivalents between
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
   combination for that index (in a boolean vector). Notice that 
   the algorithm DOES NOT calculate all possible lexicographical 
   combinations to get the requested combination.
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
