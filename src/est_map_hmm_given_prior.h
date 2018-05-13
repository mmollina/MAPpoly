#ifndef _est_EST_DIST_H
#define _est_EST_DIST_H
#include <Rcpp.h>
#include "hmm_elements.h"
#include "combinatorial.h"
#include <math.h>
#include <algorithm>
#define TOL 0
#define THRESHOLD 0.01

/*
 * note : RcppExport is an alias to `extern "C"` defined by Rcpp.
 *
 * It gives C calling convention to the rcpp_hello_world function so that
 * it can be called from .Call in R. Otherwise, the C++ compiler mangles the
 * name of the function and .Call can't find it.
 *
 * It is only useful to use RcppExport when the function is intended to be called
 * by .Call. See the thread http://thread.gmane.org/gmane.comp.lang.r.rcpp/649/focus=672
 * on Rcpp-devel for a misuse of RcppExport
 */

void setup_pre_calc_n_rec_cache(int m, int gam);

std::vector <double> alpha_ai_dot(std::vector<double>& pre_calc_prob, int m, int gam, int gam_pow_2, int l1, int l2);

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
                                 SEXP tol_R);

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
                                    SEXP verbose_R);

#endif
