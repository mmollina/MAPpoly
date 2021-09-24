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
#ifndef HMM_ELEMENTS_H
#define HMM_ELEMENTS_H

std::vector <double>  emit_poly(int m, int cte, int ip_k, int ip_k1,
                                int iq_k, int iq_k1,
                                std::vector<int>& p,
                                std::vector<int>& q,
                                std::vector<double>& g);

std::vector<std::vector<double> > transition(int m, double rf);

double prob_k1_given_k_l_m(int m, int l, double rf);

std::vector<double> forward(int m,
                            std::vector<double>& fk,
                            std::vector<int>& ik,
                            std::vector<int>& ik1,
                            std::vector<std::vector<double> >& T);

std::vector<double> backward(int m,
                             std::vector<double>& fk1,
                             std::vector<int>& ik,
                             std::vector<int>& ik1,
                             std::vector<std::vector<double> >& T);

std::vector<double> forward_emit(int m,
                                 std::vector<double>& fk,
                                 std::vector<int>& ik,
                                 std::vector<int>& ik1,
                                 std::vector<double>& emit,
                                 std::vector<std::vector<double> >& T);

std::vector<double> backward_emit(int m,
                                  std::vector<double>& fk1,
                                  std::vector<int>& ik,
                                  std::vector<int>& ik1,
                                  std::vector<double>& emit,
                                  std::vector<std::vector<double> >& T);


std::vector<double> forward_emit_one_parent(int m,
                                            std::vector<double>& fk,
                                            std::vector<int>& ik,
                                            std::vector<int>& ik1,
                                            std::vector<double>& emit,
                                            std::vector<std::vector<double> >& T);

std::vector<double> backward_emit_one_parent(int m,
                                             std::vector<double>& fk1,
                                             std::vector<int>& ik,
                                             std::vector<int>& ik1,
                                             std::vector<double>& emit,
                                             std::vector<std::vector<double> >& T);

std::vector<long double> forward_highprec(int m,
				     std::vector<long double>& fk,
				     std::vector<int>& ik,
				     std::vector<int>& ik1,
				     std::vector<std::vector<double> >& T);

std::vector<long double> backward_highprec(int m,
				      std::vector<long double>& fk1,
				      std::vector<int>& ik,
				      std::vector<int>& ik1,
				      std::vector<std::vector<double> >& T);


std::vector<long double> forward_emit_highprec(int m,
                                               std::vector<long double>& fk,
                                               std::vector<int>& ik,
                                               std::vector<int>& ik1,
                                               std::vector<double>& emit,
                                               std::vector<std::vector<double> >& T);
std::vector<long double> backward_emit_highprec(int m,
                                                std::vector<long double>& fk1,
                                                std::vector<int>& ik,
                                                std::vector<int>& ik1,
                                                std::vector<double>& emit,
                                                std::vector<std::vector<double> >& T);

std::vector<std::vector<int> > index_func(int m,
                                             std::vector<int>& p,
                                             std::vector<int>& q);

std::vector<std::vector<double> > rec_num(int m);

std::vector<std::vector<int> > rec_num_no_denominator(int m);

double init_poly(int m, int dP, int dQ, int dG);

#endif
