#include <Rcpp.h>

using namespace std;
using namespace Rcpp;

double twopt_likelihood_dosage(double rf, int m, int n_ind,
			       int dP, int dQ,
			       Rcpp::NumericVector dk,
			       Rcpp::NumericVector dk1,
			       Rcpp::NumericVector gen_1,
			       Rcpp::NumericVector gen_2,
			       Rcpp::NumericMatrix count_mat);

RcppExport SEXP pairwise_rf_estimation(SEXP m_R,
				       SEXP mrk_pairs_R,
				       SEXP geno_R,
				       SEXP dP_R,
				       SEXP dQ_R,
				       SEXP count_cache_R,
				       SEXP tol_R);
