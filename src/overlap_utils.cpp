#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>

using namespace Rcpp;

// Fast overlap detection for genomic ranges
// [[Rcpp::export]]
List fast_find_overlaps(IntegerVector start1, IntegerVector end1, 
                        IntegerVector start2, IntegerVector end2) {
    int n1 = start1.size();
    int n2 = start2.size();
    
    std::vector<int> query_hits;
    std::vector<int> subject_hits;
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            // Check if ranges overlap
            if (!(end1[i] < start2[j] || end2[j] < start1[i])) {
                query_hits.push_back(i + 1);  // R uses 1-based indexing
                subject_hits.push_back(j + 1);
            }
        }
    }
    
    return List::create(
        Named("queryHits") = wrap(query_hits),
        Named("subjectHits") = wrap(subject_hits)
    );
}

// Fast overlap detection with "within" type
// [[Rcpp::export]]
List fast_find_overlaps_within(IntegerVector start1, IntegerVector end1, 
                              IntegerVector start2, IntegerVector end2) {
    int n1 = start1.size();
    int n2 = start2.size();
    
    std::vector<int> query_hits;
    std::vector<int> subject_hits;
    
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) {
            // Check if range1 is within range2
            if (start1[i] >= start2[j] && end1[i] <= end2[j]) {
                query_hits.push_back(i + 1);
                subject_hits.push_back(j + 1);
            }
        }
    }
    
    return List::create(
        Named("queryHits") = wrap(query_hits),
        Named("subjectHits") = wrap(subject_hits)
    );
}

// Fast duplicate detection and removal for BLAST results
// [[Rcpp::export]]
DataFrame fast_remove_duplicates(DataFrame df, StringVector id_col) {
    int n = df.nrows();
    std::unordered_set<std::string> seen_ids;
    std::vector<int> keep_indices;
    
    CharacterVector ids = df[id_col[0]];
    
    for (int i = 0; i < n; i++) {
        std::string id = as<std::string>(ids[i]);
        if (seen_ids.find(id) == seen_ids.end()) {
            seen_ids.insert(id);
            keep_indices.push_back(i);
        }
    }
    
    // Create new dataframe with kept rows
    DataFrame result = DataFrame::create();
    for (int j = 0; j < df.size(); j++) {
        result.push_back(df[j]);
    }
    
    // Filter rows
    for (int j = 0; j < result.size(); j++) {
        SEXP col = result[j];
        if (TYPEOF(col) == INTSXP) {
            IntegerVector new_col(keep_indices.size());
            for (size_t k = 0; k < keep_indices.size(); k++) {
                new_col[k] = as<IntegerVector>(col)[keep_indices[k]];
            }
            result[j] = new_col;
        } else if (TYPEOF(col) == REALSXP) {
            NumericVector new_col(keep_indices.size());
            for (size_t k = 0; k < keep_indices.size(); k++) {
                new_col[k] = as<NumericVector>(col)[keep_indices[k]];
            }
            result[j] = new_col;
        } else if (TYPEOF(col) == STRSXP) {
            CharacterVector new_col(keep_indices.size());
            for (size_t k = 0; k < keep_indices.size(); k++) {
                new_col[k] = as<CharacterVector>(col)[keep_indices[k]];
            }
            result[j] = new_col;
        }
    }
    
    return result;
}

// Fast sorting and grouping for RBH operations
// [[Rcpp::export]]
DataFrame fast_sort_by_column(DataFrame df, StringVector sort_col, bool decreasing = true) {
    int n = df.nrows();
    std::vector<std::pair<double, int>> sort_pairs;
    
    NumericVector sort_values = df[sort_col[0]];
    
    for (int i = 0; i < n; i++) {
        sort_pairs.push_back(std::make_pair(sort_values[i], i));
    }
    
    if (decreasing) {
        std::sort(sort_pairs.begin(), sort_pairs.end(), 
                 std::greater<std::pair<double, int>>());
    } else {
        std::sort(sort_pairs.begin(), sort_pairs.end());
    }
    
    // Create new dataframe with sorted rows
    DataFrame result = DataFrame::create();
    for (int j = 0; j < df.size(); j++) {
        result.push_back(df[j]);
    }
    
    // Reorder rows
    for (int j = 0; j < result.size(); j++) {
        SEXP col = result[j];
        if (TYPEOF(col) == INTSXP) {
            IntegerVector new_col(n);
            for (int k = 0; k < n; k++) {
                new_col[k] = as<IntegerVector>(col)[sort_pairs[k].second];
            }
            result[j] = new_col;
        } else if (TYPEOF(col) == REALSXP) {
            NumericVector new_col(n);
            for (int k = 0; k < n; k++) {
                new_col[k] = as<NumericVector>(col)[sort_pairs[k].second];
            }
            result[j] = new_col;
        } else if (TYPEOF(col) == STRSXP) {
            CharacterVector new_col(n);
            for (int k = 0; k < n; k++) {
                new_col[k] = as<CharacterVector>(col)[sort_pairs[k].second];
            }
            result[j] = new_col;
        }
    }
    
    return result;
}
