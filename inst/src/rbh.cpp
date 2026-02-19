#include <Rcpp.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

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
};

// parse 8 columns from BLAST outfmt6 line:
// query_tx subject_tx pident qcovs qstart qend sstart send
static inline bool parse_blast8(const std::string &line,
                                std::string &q, std::string &s,
                                double &pident, int &qcovs,
                                int &qstart, int &qend,
                                int &sstart, int &send){
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
            else if(col==8) { send = std::atoi(p); return true; }
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
// txA txB Astart Aend Bstart Bend pident dir qcovs
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
    int qcovs=0, qstart=0, qend=0, sstart=0, send=0;
    
    while(std::getline(in, line)){
        trim_cr(line);
        if(!parse_blast8(line, q, s, pident, qcovs, qstart, qend, sstart, send)) continue;
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
        
        if(Astart > Aend) std::swap(Astart, Aend);
        if(Bstart > Bend) std::swap(Bstart, Bend);
        
        out << txA << '\t' << txB << '\t'
            << Astart << '\t' << Aend << '\t'
            << Bstart << '\t' << Bend << '\t'
            << pident << '\t' << dir << '\t' << qcovs
            << '\n';
    }
    out.close(); in.close();
}

static bool read_norm_line(std::ifstream &in, NormRec &r){
    std::string line;
    if(!std::getline(in, line)) return false;
    trim_cr(line);
    
    // 9 columns: txA txB Astart Aend Bstart Bend pident dir qcovs
    size_t tab[8];
    size_t pos=0;
    for(int i=0;i<8;i++){
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
    return true;
}

// Stream sorted normalized file and write RBH table.
// Reciprocal criterion: (txA, txB) has at least one HSP with dir=1 and one with dir=2.
// pident = length-weighted mean over ALL HSPs for the pair.
// q2s_qcovs = max qcovs over dir=1, s2q_qcovs = max qcovs over dir=2.
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
    bool has_dir1 = false, has_dir2 = false;
    int q2s_qcovs = NA_INTEGER, s2q_qcovs = NA_INTEGER;
    double ident_sum = 0.0;
    long long aln_sum = 0;
    
    auto flush_pair = [&](){
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
        if(cur.dir == 1){
            has_dir1 = true;
            if(q2s_qcovs == NA_INTEGER || cur.qcovs > q2s_qcovs) q2s_qcovs = cur.qcovs;
        } else if(cur.dir == 2){
            has_dir2 = true;
            if(s2q_qcovs == NA_INTEGER || cur.qcovs > s2q_qcovs) s2q_qcovs = cur.qcovs;
        }
        int alen = cur.Aend - cur.Astart + 1;
        if(alen > 0){
            aln_sum += (long long)alen;
            ident_sum += (cur.pident * 0.01) * (double)alen;
        }
        
        if(!read_norm_line(in, nxt)){
            flush_pair();
            out.close(); in.close();
            return;
        }
        
        if(nxt.txA != pairA || nxt.txB != pairB){
            flush_pair();
            pairA = nxt.txA; pairB = nxt.txB;
            has_dir1 = false; has_dir2 = false;
            q2s_qcovs = NA_INTEGER; s2q_qcovs = NA_INTEGER;
            ident_sum = 0.0; aln_sum = 0;
        }
        cur = nxt;
    }
}
