#ifndef COMBINATORIAL_FUNCTIONS_H
#define COMBINATORIAL_FUNCTIONS_H

int nChoosek(int n, int k);
int n_rec_given_genk_and_k1(int ploidy, int index1, int index2);
double prob_k1_given_k_lp_lq_m(int m, int lp, int lq, double rf);
std::vector <int>  boolean_lexicographic_k_choose_m_and_collapse(int ploidy,
        std::vector<int>& which_homologous_mk1,
        std::vector<int>& which_homologous_mk2,
        int gen_prog_mk1,
        int gen_prog_mk2);
std::vector <bool> get_boolean_vec_from_lexicographical_index(int ploidy, int index);
#endif

