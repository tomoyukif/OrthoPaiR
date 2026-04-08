#include <Rcpp.h>
#include <unordered_map>
#include <vector>
#include <string>
using namespace Rcpp;

// [[Rcpp::export]]
List rbh_shortest_path(
        IntegerVector qgeneid,
        IntegerVector sgeneid,
        IntegerVector q_chr,
        IntegerVector s_chr
) {
    int n = qgeneid.size();
    
    if (sgeneid.size() != n || q_chr.size() != n || s_chr.size() != n) {
        stop("All input vectors must have the same length.");
    }
    
    IntegerVector out_idx(n, NA_INTEGER);
    NumericVector out_dist(n, NA_REAL);
    
    std::unordered_map<std::string, std::vector<int>> groups;
    groups.reserve(n);
    
    // グループ化
    for (int i = 0; i < n; ++i) {
        std::string key = std::to_string(q_chr[i]) + "|" + std::to_string(s_chr[i]);
        groups[key].push_back(i);
    }
    
    // 各グループで処理
    for (auto &kv : groups) {
        std::vector<int> &idx = kv.second;
        int m = idx.size();
        
        if (m <= 1) {
            out_idx[idx[0]] = NA_INTEGER;
            out_dist[idx[0]] = NA_REAL;
            continue;
        }
        
        for (int a = 0; a < m; ++a) {
            int i = idx[a];
            
            long long best_dist = LLONG_MAX;
            int best_idx = NA_INTEGER;
            
            int qi = qgeneid[i];
            int si = sgeneid[i];
            
            for (int b = 0; b < m; ++b) {
                int j = idx[b];
                if (i == j) continue;
                
                long long dq = (long long)qgeneid[j] - qi;
                long long ds = (long long)sgeneid[j] - si;
                
                // 元コード準拠: qgeneid同一 or sgeneid同一 は無効
                if (dq == 0 || ds == 0) continue;
                
                long long dist = std::llabs(dq) + std::llabs(ds);
                
                if (dist < best_dist) {
                    best_dist = dist;
                    best_idx = j + 1;   // R index
                }
            }
            
            out_idx[i] = best_idx;
            out_dist[i] = (double) best_dist;
        }
    }
    
    return List::create(
        _["shortest_path"] = out_idx,
        _["best_dist"] = out_dist
    );
}