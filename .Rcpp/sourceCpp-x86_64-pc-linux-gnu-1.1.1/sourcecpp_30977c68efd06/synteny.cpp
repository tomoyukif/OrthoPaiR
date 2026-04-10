#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
#include <climits>
using namespace Rcpp;

// For each RBH row, find the nearest anchor gene pair (same query/subject chromosome
// pair) minimizing |anchor_q - q| + |anchor_s - s|. Returns 1-based index into the
// anchor table and the L1 distance.
// If omit_zero is true, anchor pairs at distance 0 are ignored (returns NA if none left).

// [[Rcpp::export]]
List rbh_shortest_path(
        IntegerVector qgeneid,
        IntegerVector sgeneid,
        IntegerVector q_chr,
        IntegerVector s_chr,
        IntegerVector anchor_qgeneid,
        IntegerVector anchor_sgeneid,
        IntegerVector anchor_q_chr,
        IntegerVector anchor_s_chr,
        bool omit_zero = false
) {
    int n = qgeneid.size();
    
    if (sgeneid.size() != n || q_chr.size() != n || s_chr.size() != n) {
        stop("RBH vectors (qgeneid, sgeneid, q_chr, s_chr) must have the same length.");
    }
    
    int na = anchor_qgeneid.size();
    if (na == 0) {
        stop("anchor_* vectors must be non-empty.");
    }
    if ((int) anchor_sgeneid.size() != na ||
        (int) anchor_q_chr.size() != na ||
        (int) anchor_s_chr.size() != na) {
        stop("All anchor_* vectors must have the same length.");
    }
    
    IntegerVector out_idx(n, NA_INTEGER);
    NumericVector out_dist(n, NA_REAL);
    
    std::unordered_map<std::string, std::vector<int>> anchor_by_chrpair;
    anchor_by_chrpair.reserve((size_t) na);
    
    for (int k = 0; k < na; ++k) {
        std::string key = std::to_string(anchor_q_chr[k]) + "|" + std::to_string(anchor_s_chr[k]);
        anchor_by_chrpair[key].push_back(k);
    }
    
    for (int i = 0; i < n; ++i) {
        std::string key = std::to_string(q_chr[i]) + "|" + std::to_string(s_chr[i]);
        auto it = anchor_by_chrpair.find(key);
        if (it == anchor_by_chrpair.end() || it->second.empty()) {
            out_idx[i] = NA_INTEGER;
            out_dist[i] = NA_REAL;
            continue;
        }
        
        const std::vector<int> &cand = it->second;
        long long best_dist = LLONG_MAX;
        int best_anchor_1based = NA_INTEGER;
        
        int qi = qgeneid[i];
        int si = sgeneid[i];
        
        for (int t = 0; t < (int) cand.size(); ++t) {
            int k = cand[t];
            long long dq = (long long) anchor_qgeneid[k] - qi;
            long long ds = (long long) anchor_sgeneid[k] - si;
            long long dist = std::llabs(dq) + std::llabs(ds);
            
            if (omit_zero && dist == 0) {
                continue;
            }
            if (dist < best_dist) {
                best_dist = dist;
                best_anchor_1based = k + 1;  // R index into anchor table
            }
        }
        
        if (best_dist == LLONG_MAX) {
            out_idx[i] = NA_INTEGER;
            out_dist[i] = NA_REAL;
        } else {
            out_idx[i] = best_anchor_1based;
            out_dist[i] = (double) best_dist;
        }
    }
    
    return List::create(
        _["shortest_path"] = out_idx,
        _["best_dist"] = out_dist
    );
}


#include <Rcpp.h>
#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rbh_shortest_path
List rbh_shortest_path(IntegerVector qgeneid, IntegerVector sgeneid, IntegerVector q_chr, IntegerVector s_chr, IntegerVector anchor_qgeneid, IntegerVector anchor_sgeneid, IntegerVector anchor_q_chr, IntegerVector anchor_s_chr, bool omit_zero);
RcppExport SEXP sourceCpp_1_rbh_shortest_path(SEXP qgeneidSEXP, SEXP sgeneidSEXP, SEXP q_chrSEXP, SEXP s_chrSEXP, SEXP anchor_qgeneidSEXP, SEXP anchor_sgeneidSEXP, SEXP anchor_q_chrSEXP, SEXP anchor_s_chrSEXP, SEXP omit_zeroSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type qgeneid(qgeneidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type sgeneid(sgeneidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type q_chr(q_chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type s_chr(s_chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anchor_qgeneid(anchor_qgeneidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anchor_sgeneid(anchor_sgeneidSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anchor_q_chr(anchor_q_chrSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type anchor_s_chr(anchor_s_chrSEXP);
    Rcpp::traits::input_parameter< bool >::type omit_zero(omit_zeroSEXP);
    rcpp_result_gen = Rcpp::wrap(rbh_shortest_path(qgeneid, sgeneid, q_chr, s_chr, anchor_qgeneid, anchor_sgeneid, anchor_q_chr, anchor_s_chr, omit_zero));
    return rcpp_result_gen;
END_RCPP
}
