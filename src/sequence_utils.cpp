#include <Rcpp.h>
#include <vector>
#include <algorithm>
#include <string>
#include <unordered_map>

using namespace Rcpp;

// Fast sequence similarity calculation using simple scoring
// [[Rcpp::export]]
double fast_sequence_similarity(std::string seq1, std::string seq2) {
    int len1 = seq1.length();
    int len2 = seq2.length();
    
    if (len1 == 0 || len2 == 0) return 0.0;
    
    // Simple alignment scoring (can be enhanced)
    int matches = 0;
    int min_len = std::min(len1, len2);
    
    for (int i = 0; i < min_len; i++) {
        if (seq1[i] == seq2[i]) {
            matches++;
        }
    }
    
    return (double)matches / std::max(len1, len2);
}

// Fast pairwise alignment scoring for multiple sequences
// [[Rcpp::export]]
NumericMatrix fast_pairwise_scores(CharacterVector sequences) {
    int n = sequences.size();
    NumericMatrix scores(n, n);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (i == j) {
                scores(i, j) = 1.0;
            } else {
                std::string seq1 = as<std::string>(sequences[i]);
                std::string seq2 = as<std::string>(sequences[j]);
                scores(i, j) = fast_sequence_similarity(seq1, seq2);
            }
        }
    }
    
    return scores;
}

// Fast grouping of similar sequences
// [[Rcpp::export]]
List fast_group_similar_sequences(CharacterVector sequences, double threshold = 0.8) {
    int n = sequences.size();
    std::vector<int> group_ids(n, -1);
    int current_group = 0;
    
    for (int i = 0; i < n; i++) {
        if (group_ids[i] == -1) {
            group_ids[i] = current_group;
            
            for (int j = i + 1; j < n; j++) {
                if (group_ids[j] == -1) {
                    std::string seq1 = as<std::string>(sequences[i]);
                    std::string seq2 = as<std::string>(sequences[j]);
                    double similarity = fast_sequence_similarity(seq1, seq2);
                    
                    if (similarity >= threshold) {
                        group_ids[j] = current_group;
                    }
                }
            }
            current_group++;
        }
    }
    
    return List::create(
        Named("group_ids") = wrap(group_ids),
        Named("n_groups") = current_group
    );
}

// Fast string operations for ID manipulation
// [[Rcpp::export]]
CharacterVector fast_paste_strings(CharacterVector vec1, CharacterVector vec2, 
                                  std::string separator = "_") {
    int n = vec1.size();
    CharacterVector result(n);
    
    for (int i = 0; i < n; i++) {
        std::string str1 = as<std::string>(vec1[i]);
        std::string str2 = as<std::string>(vec2[i]);
        result[i] = str1 + separator + str2;
    }
    
    return result;
}

// Fast ID replacement using hash map
// [[Rcpp::export]]
CharacterVector fast_replace_ids(CharacterVector ids, CharacterVector old_ids, 
                               CharacterVector new_ids) {
    int n = ids.size();
    int n_replace = old_ids.size();
    
    // Create hash map for fast lookup
    std::unordered_map<std::string, std::string> id_map;
    for (int i = 0; i < n_replace; i++) {
        id_map[as<std::string>(old_ids[i])] = as<std::string>(new_ids[i]);
    }
    
    CharacterVector result(n);
    for (int i = 0; i < n; i++) {
        std::string id = as<std::string>(ids[i]);
        auto it = id_map.find(id);
        if (it != id_map.end()) {
            result[i] = it->second;
        } else {
            result[i] = id;
        }
    }
    
    return result;
}

// Fast counting of unique values
// [[Rcpp::export]]
List fast_count_unique(CharacterVector vec) {
    std::unordered_map<std::string, int> counts;
    int n = vec.size();
    
    for (int i = 0; i < n; i++) {
        std::string val = as<std::string>(vec[i]);
        counts[val]++;
    }
    
    CharacterVector unique_vals;
    IntegerVector counts_vec;
    
    for (auto& pair : counts) {
        unique_vals.push_back(pair.first);
        counts_vec.push_back(pair.second);
    }
    
    return List::create(
        Named("values") = unique_vals,
        Named("counts") = counts_vec
    );
}
