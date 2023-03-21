// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// pairwise_rf_estimation_disc_rcpp
RcppExport SEXP pairwise_rf_estimation_disc_rcpp(SEXP mrk_pairs_R, SEXP m_R, SEXP geno_R, SEXP dP_R, SEXP dQ_R, SEXP count_vector_R, SEXP count_matrix_phases_R, SEXP count_matrix_rownames_R, SEXP count_matrix_number_R, SEXP count_matrix_pos_R, SEXP count_matrix_length_R, SEXP tol_R, SEXP threads_R);
RcppExport SEXP _mappoly_pairwise_rf_estimation_disc_rcpp(SEXP mrk_pairs_RSEXP, SEXP m_RSEXP, SEXP geno_RSEXP, SEXP dP_RSEXP, SEXP dQ_RSEXP, SEXP count_vector_RSEXP, SEXP count_matrix_phases_RSEXP, SEXP count_matrix_rownames_RSEXP, SEXP count_matrix_number_RSEXP, SEXP count_matrix_pos_RSEXP, SEXP count_matrix_length_RSEXP, SEXP tol_RSEXP, SEXP threads_RSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< SEXP >::type mrk_pairs_R(mrk_pairs_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type m_R(m_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type geno_R(geno_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dP_R(dP_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type dQ_R(dQ_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_vector_R(count_vector_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_matrix_phases_R(count_matrix_phases_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_matrix_rownames_R(count_matrix_rownames_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_matrix_number_R(count_matrix_number_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_matrix_pos_R(count_matrix_pos_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type count_matrix_length_R(count_matrix_length_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type tol_R(tol_RSEXP);
    Rcpp::traits::input_parameter< SEXP >::type threads_R(threads_RSEXP);
    rcpp_result_gen = Rcpp::wrap(pairwise_rf_estimation_disc_rcpp(mrk_pairs_R, m_R, geno_R, dP_R, dQ_R, count_vector_R, count_matrix_phases_R, count_matrix_rownames_R, count_matrix_number_R, count_matrix_pos_R, count_matrix_length_R, tol_R, threads_R));
    return rcpp_result_gen;
END_RCPP
}
// vcf_get_probabilities
Rcpp::List vcf_get_probabilities(Rcpp::StringMatrix& mat, int pl_pos);
RcppExport SEXP _mappoly_vcf_get_probabilities(SEXP matSEXP, SEXP pl_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type pl_pos(pl_posSEXP);
    rcpp_result_gen = Rcpp::wrap(vcf_get_probabilities(mat, pl_pos));
    return rcpp_result_gen;
END_RCPP
}
// vcf_transform_dosage
Rcpp::NumericMatrix vcf_transform_dosage(Rcpp::StringMatrix& mat, int gt_pos);
RcppExport SEXP _mappoly_vcf_transform_dosage(SEXP matSEXP, SEXP gt_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type gt_pos(gt_posSEXP);
    rcpp_result_gen = Rcpp::wrap(vcf_transform_dosage(mat, gt_pos));
    return rcpp_result_gen;
END_RCPP
}
// vcf_get_ploidy
Rcpp::NumericMatrix vcf_get_ploidy(Rcpp::StringMatrix& mat, int gt_pos);
RcppExport SEXP _mappoly_vcf_get_ploidy(SEXP matSEXP, SEXP gt_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type gt_pos(gt_posSEXP);
    rcpp_result_gen = Rcpp::wrap(vcf_get_ploidy(mat, gt_pos));
    return rcpp_result_gen;
END_RCPP
}
// vcf_get_depth
Rcpp::NumericMatrix vcf_get_depth(Rcpp::StringMatrix& mat, int dp_pos);
RcppExport SEXP _mappoly_vcf_get_depth(SEXP matSEXP, SEXP dp_posSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::StringMatrix& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< int >::type dp_pos(dp_posSEXP);
    rcpp_result_gen = Rcpp::wrap(vcf_get_depth(mat, dp_pos));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP calc_genoprob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP calc_genoprob_prior(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP calc_genprob_haplo(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP calc_genprob_haplo_highprec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP calc_genprob_single_parent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_haplotype_map(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_haplotype_map_highprec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_hmm_map_single_parent(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_map_hmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP est_map_hmm_highprec(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP get_counts_single_parent_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP loglike_hmm(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP pairwise_rf_estimation_disc(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP pairwise_rf_estimation_disc_rcpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP pairwise_rf_estimation_prob(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
RcppExport SEXP poly_hmm_est_CPP(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_mappoly_pairwise_rf_estimation_disc_rcpp", (DL_FUNC) &_mappoly_pairwise_rf_estimation_disc_rcpp, 13},
    {"_mappoly_vcf_get_probabilities", (DL_FUNC) &_mappoly_vcf_get_probabilities, 2},
    {"_mappoly_vcf_transform_dosage", (DL_FUNC) &_mappoly_vcf_transform_dosage, 2},
    {"_mappoly_vcf_get_ploidy", (DL_FUNC) &_mappoly_vcf_get_ploidy, 2},
    {"_mappoly_vcf_get_depth", (DL_FUNC) &_mappoly_vcf_get_depth, 2},
    {"calc_genoprob",                    (DL_FUNC) &calc_genoprob,                     7},
    {"calc_genoprob_prior",              (DL_FUNC) &calc_genoprob_prior,              12},
    {"calc_genprob_haplo",               (DL_FUNC) &calc_genprob_haplo,                8},
    {"calc_genprob_haplo_highprec",      (DL_FUNC) &calc_genprob_haplo_highprec,       8},
    {"calc_genprob_single_parent",       (DL_FUNC) &calc_genprob_single_parent,        8},
    {"est_haplotype_map",                (DL_FUNC) &est_haplotype_map,                 9},
    {"est_haplotype_map_highprec",       (DL_FUNC) &est_haplotype_map_highprec,        9},
    {"est_hmm_map_single_parent",        (DL_FUNC) &est_hmm_map_single_parent,         9},
    {"est_map_hmm",                      (DL_FUNC) &est_map_hmm,                       8},
    {"est_map_hmm_highprec",             (DL_FUNC) &est_map_hmm_highprec,              8},
    {"get_counts_single_parent_cpp",     (DL_FUNC) &get_counts_single_parent_cpp,      6},
    {"loglike_hmm",                      (DL_FUNC) &loglike_hmm,                       6},
    {"pairwise_rf_estimation_disc",      (DL_FUNC) &pairwise_rf_estimation_disc,       7},
    {"pairwise_rf_estimation_disc_rcpp", (DL_FUNC) &pairwise_rf_estimation_disc_rcpp, 13},
    {"pairwise_rf_estimation_prob",      (DL_FUNC) &pairwise_rf_estimation_prob,       8},
    {"poly_hmm_est_CPP",                 (DL_FUNC) &poly_hmm_est_CPP,                 13},
    {NULL, NULL, 0}
};

RcppExport void R_init_mappoly(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
