#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <set>
#include <utility>

using namespace Rcpp;

static inline void trim_cr(std::string &s){
    if(!s.empty() && s.back()=='\r') s.pop_back();
}

struct NormRec {
    std::string txA, txB;
    int Astart, Aend, Bstart, Bend;
    double pident; // keep but NOT part of reciprocal key anymore
    int dir;       // 1 or 2
    int qcovs;
    int len_A;     // transcript length of txA (from BLAST qlen/slen)
    int len_B;     // transcript length of txB
};

// Parse 10 columns from BLAST outfmt6 line:
// qseqid sseqid pident qcovs qstart qend sstart send qlen slen
static inline bool parse_blast10(const std::string &line,
                                 std::string &q, std::string &s,
                                 double &pident, int &qcovs,
                                 int &qstart, int &qend,
                                 int &sstart, int &send,
                                 int &qlen, int &slen){
    int col = 1;
    size_t start = 0;
    size_t n = line.size();
    for(size_t i=0;i<=n;i++){
        if(i==n || line[i]=='\t'){
            const char* p = line.c_str() + start;
            if(col==1) q.assign(line, start, i-start);
            else if(col==2) s.assign(line, start, i-start);
            else if(col==3) pident = std::strtod(p, nullptr);
            else if(col==4) qcovs  = std::atoi(p);
            else if(col==5) qstart = std::atoi(p);
            else if(col==6) qend   = std::atoi(p);
            else if(col==7) sstart = std::atoi(p);
            else if(col==8) send   = std::atoi(p);
            else if(col==9) qlen   = std::atoi(p);
            else if(col==10) { slen = std::atoi(p); return true; }
            col++;
            start = i+1;
        }
    }
    return false;
}

static inline int genome_prefix_int(const std::string &tx, int width){
    if((int)tx.size() < width) return 0;
    return std::atoi(tx.substr(0, width).c_str());
}

// Normalized TSV columns:
// txA txB Astart Aend Bstart Bend pident dir qcovs len_A len_B
// [[Rcpp::export]]
void normalize_blast_tsv_cpp(std::string infile,
                             std::string outfile,
                             int genome_width = 4,
                             bool drop_self = true){
    std::ifstream in(infile, std::ios::binary);
    if(!in) stop("Cannot open infile: %s", infile);
    std::ofstream out(outfile, std::ios::binary);
    if(!out) stop("Cannot open outfile: %s", outfile);
    
    std::string line, q, s;
    double pident=0.0;
    int qcovs=0, qstart=0, qend=0, sstart=0, send=0, qlen=0, slen=0;
    
    while(std::getline(in, line)){
        trim_cr(line);
        if(!parse_blast10(line, q, s, pident, qcovs, qstart, qend, sstart, send, qlen, slen)) continue;
        if(drop_self && q==s) continue;
        
        int qg = genome_prefix_int(q, genome_width);
        int sg = genome_prefix_int(s, genome_width);
        
        // For inter-genome hits: normalize by genome (smaller genome first)
        // For intra-genome hits (qg == sg): normalize by transcript ID to preserve direction
        bool flip;
        if(qg != sg){
            flip = (qg > sg);  // Inter-genome: flip if query genome > subject genome
        } else {
            flip = (q > s);    // Intra-genome: flip if query transcript > subject transcript
        }
        std::string txA = flip ? s : q;
        std::string txB = flip ? q : s;
        int dir = flip ? 2 : 1;
        
        int Astart = flip ? sstart : qstart;
        int Aend   = flip ? send   : qend;
        int Bstart = flip ? qstart : sstart;
        int Bend   = flip ? qend   : send;
        // len_A = length of txA, len_B = length of txB (unchanged by flip)
        int len_A = flip ? slen : qlen;
        int len_B = flip ? qlen : slen;
        
        if(Astart > Aend) std::swap(Astart, Aend);
        if(Bstart > Bend) std::swap(Bstart, Bend);
        
        out << txA << '\t' << txB << '\t'
            << Astart << '\t' << Aend << '\t'
            << Bstart << '\t' << Bend << '\t'
            << pident << '\t' << dir << '\t' << qcovs << '\t'
            << len_A << '\t' << len_B << '\n';
    }
    out.close(); in.close();
}

static bool read_norm_line(std::ifstream &in, NormRec &r){
    std::string line;
    if(!std::getline(in, line)) return false;
    trim_cr(line);
    
    // 11 columns: txA txB Astart Aend Bstart Bend pident dir qcovs len_A len_B
    size_t tab[10];
    size_t pos=0;
    for(int i=0;i<10;i++){
        size_t t = line.find('\t', pos);
        if(t==std::string::npos) return false;
        tab[i]=t;
        pos=t+1;
    }
    r.txA = line.substr(0, tab[0]);
    r.txB = line.substr(tab[0]+1, tab[1]-tab[0]-1);
    
    r.Astart = std::atoi(line.c_str()+tab[1]+1);
    r.Aend   = std::atoi(line.c_str()+tab[2]+1);
    r.Bstart = std::atoi(line.c_str()+tab[3]+1);
    r.Bend   = std::atoi(line.c_str()+tab[4]+1);
    r.pident = std::strtod(line.c_str()+tab[5]+1, nullptr);
    r.dir    = std::atoi(line.c_str()+tab[6]+1);
    r.qcovs  = std::atoi(line.c_str()+tab[7]+1);
    r.len_A  = std::atoi(line.c_str()+tab[8]+1);
    r.len_B  = std::atoi(line.c_str()+tab[9]+1);
    return true;
}

// r1 is contained in r2 on query (A) or subject (B) — can be equal (duplicate)
static inline bool is_contained_or_equal_query(const NormRec &r1, const NormRec &r2){
    return (r2.Astart <= r1.Astart && r1.Aend <= r2.Aend);
}
static inline bool is_contained_or_equal_subject(const NormRec &r1, const NormRec &r2){
    return (r2.Bstart <= r1.Bstart && r1.Bend <= r2.Bend);
}

// r1 and r2 have identical alignment ranges (A and B)
static inline bool ranges_equal(const NormRec &r1, const NormRec &r2){
    return (r1.Astart == r2.Astart && r1.Aend == r2.Aend &&
            r1.Bstart == r2.Bstart && r1.Bend == r2.Bend);
}

// Filter: remove HSPs that are contained or duplicate only within the same direction.
// - Different dir (1 vs 2): bidirectional hit → keep both.
// - Same dir: remove HSPs that are contained in another, or duplicate (identical range);
//   keep one representative per group (first occurrence).
static std::vector<size_t> filter_contained_hsps(const std::vector<NormRec> &hsps){
    std::vector<bool> is_contained(hsps.size(), false);
    
    for(size_t i = 0; i < hsps.size(); i++){
        for(size_t j = 0; j < hsps.size(); j++){
            if(i == j) continue;
            if(hsps[i].dir != hsps[j].dir) continue; // different dir = bidirectional, keep both
            // Same dir: remove i if (1) i is strictly contained in j, or (2) duplicate (keep first)
            bool duplicate = ranges_equal(hsps[i], hsps[j]) && (i > j);
            bool contained = !ranges_equal(hsps[i], hsps[j]) &&
                (is_contained_or_equal_query(hsps[i], hsps[j]) ||
                 is_contained_or_equal_subject(hsps[i], hsps[j]));
            if(contained || duplicate){
                is_contained[i] = true;
                break;
            }
        }
    }
    
    std::vector<size_t> keep;
    keep.reserve(hsps.size());
    for(size_t i = 0; i < hsps.size(); i++){
        if(!is_contained[i]) keep.push_back(i);
    }
    return keep;
}

// Calculate union length of alignment ranges
// Ranges are represented as (start, end) pairs
static int union_length(const std::vector<std::pair<int, int> > &ranges){
    if(ranges.empty()) return 0;
    
    // Sort ranges by start position
    std::vector<std::pair<int, int> > sorted = ranges;
    std::sort(sorted.begin(), sorted.end());
    
    int total = 0;
    int current_start = sorted[0].first;
    int current_end = sorted[0].second;
    
    for(size_t i = 1; i < sorted.size(); i++){
        if(sorted[i].first <= current_end){
            // Overlapping or adjacent: extend current range
            if(sorted[i].second > current_end) current_end = sorted[i].second;
        } else {
            // Non-overlapping: add current range and start new one
            total += current_end - current_start + 1;
            current_start = sorted[i].first;
            current_end = sorted[i].second;
        }
    }
    total += current_end - current_start + 1;
    
    return total;
}

// Recalculate coverage from union of alignment ranges using true transcript length.
// tx_len must come from BLAST qlen/slen (output by blastn), not from qcovs (qcovs is
// overall coverage over all HSPs for that query, not per-HSP).
static int recalc_coverage(const std::vector<NormRec> &hsps, 
                           const std::vector<size_t> &keep_idx,
                           int dir, bool use_query, int tx_len){
    if(tx_len <= 0) return NA_INTEGER;
    
    std::vector<std::pair<int, int> > ranges;
    for(size_t idx : keep_idx){
        const NormRec &r = hsps[idx];
        if(r.dir != dir) continue;
        
        int start, end;
        if(use_query){
            start = r.Astart;
            end = r.Aend;
        } else {
            start = r.Bstart;
            end = r.Bend;
        }
        ranges.push_back(std::make_pair(start, end));
    }
    
    if(ranges.empty()) return NA_INTEGER;
    
    int union_len = union_length(ranges);
    // coverage = 100 * union_len / tx_len; +0.5 for rounding to nearest integer
    return (int)((double)union_len * 100.0 / (double)tx_len + 0.5);
}

// Stream sorted normalized file and write RBH table.
// Reciprocal criterion: (txA, txB) has at least one HSP with dir=1 and one with dir=2.
// pident = length-weighted mean over ALL HSPs for the pair (after filtering contained HSPs).
// q2s_qcovs = recalculated coverage on txA from union of dir=1 HSP alignment ranges.
// s2q_qcovs = recalculated coverage on txB from union of dir=2 HSP alignment ranges.
// Output: query_tx, subject_tx, pident, q2s_qcovs, s2q_qcovs, q2s_ci, s2q_ci, mutual_ci (8 columns).
// Assumes input is sorted by txA, txB only.
// [[Rcpp::export]]
void rbh_from_sorted_norm_cpp(std::string sorted_norm_fn,
                              std::string out_fn){
    std::ifstream in(sorted_norm_fn, std::ios::binary);
    if(!in) stop("Cannot open sorted_norm_fn: %s", sorted_norm_fn);
    std::ofstream out(out_fn, std::ios::binary);
    if(!out) stop("Cannot open out_fn: %s", out_fn);
    
    NormRec cur, nxt;
    if(!read_norm_line(in, cur)){
        out.close(); in.close();
        return;
    }
    
    std::string pairA = cur.txA, pairB = cur.txB;
    std::vector<NormRec> pair_hsps;
    
    auto flush_pair = [&](){
        if(pair_hsps.empty()) return;
        
        // Filter out HSPs completely contained within another HSP
        std::vector<size_t> keep_idx = filter_contained_hsps(pair_hsps);
        if(keep_idx.empty()) return;
        
        bool has_dir1 = false, has_dir2 = false;
        double ident_sum = 0.0;
        long long aln_sum = 0;
        
        // Calculate pident from filtered HSPs
        for(size_t idx : keep_idx){
            const NormRec &r = pair_hsps[idx];
            if(r.dir == 1) has_dir1 = true;
            else if(r.dir == 2) has_dir2 = true;
            
            int alen = r.Aend - r.Astart + 1;
            if(alen > 0){
                aln_sum += (long long)alen;
                ident_sum += (r.pident * 0.01) * (double)alen;
            }
        }
        
        // Transcript lengths from BLAST qlen/slen (same for all HSPs of this pair)
        int len_A = pair_hsps[keep_idx[0]].len_A;
        int len_B = pair_hsps[keep_idx[0]].len_B;
        // Recalculate coverage from union of ranges using true transcript length
        int q2s_qcovs = recalc_coverage(pair_hsps, keep_idx, 1, true, len_A);   // coverage on txA
        int s2q_qcovs = recalc_coverage(pair_hsps, keep_idx, 2, false, len_B);   // coverage on txB
        
        if(!has_dir1 || !has_dir2 || aln_sum <= 0) return;
        double pident_pair = ident_sum / (double)aln_sum * 100.0;
        double ci_q2s = NA_REAL, ci_s2q = NA_REAL, mutual = NA_REAL;
        if(q2s_qcovs != NA_INTEGER) ci_q2s = pident_pair * (double)q2s_qcovs * 1e-4;
        if(s2q_qcovs != NA_INTEGER) ci_s2q = pident_pair * (double)s2q_qcovs * 1e-4;
        if(ci_q2s == ci_q2s && ci_s2q == ci_s2q) mutual = ci_q2s * ci_s2q;
        out << pairA << '\t' << pairB << '\t'
            << pident_pair << '\t'
            << q2s_qcovs << '\t' << s2q_qcovs << '\t'
            << ci_q2s << '\t' << ci_s2q << '\t' << mutual << '\n';
    };
    
    while(true){
        pair_hsps.push_back(cur);
        
        if(!read_norm_line(in, nxt)){
            flush_pair();
            out.close(); in.close();
            return;
        }
        
        if(nxt.txA != pairA || nxt.txB != pairB){
            flush_pair();
            pairA = nxt.txA; pairB = nxt.txB;
            pair_hsps.clear();
        }
        cur = nxt;
    }
}
