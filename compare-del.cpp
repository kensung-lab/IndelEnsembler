#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <climits>
#include <stdio.h>
#include <unistd.h>
#include "IntervalTree.h"
#include "ssw_cpp.h"
#include "kseq.h"
KSEQ_INIT(int, read)

using namespace std;
using namespace StripedSmithWaterman;

struct bp_t {
    string chr;
    int pos;

    bp_t(string chr, int pos) : chr(chr), pos(pos) {}
};

struct sv_t {
    string id;
    bp_t bp1, bp2;

    sv_t(string id, bp_t bp1, bp_t bp2) : id(id), bp1(bp1), bp2(bp2) {}

    int size() {return bp2.pos-bp1.pos;}
};

struct repeat_t {
    string chr;
    int start, end;
    
    repeat_t() {}
    repeat_t(string chr, int start, int end) : chr(chr), start(start), end(end) {}
};
bool operator == (const repeat_t& r1, const repeat_t& r2) {
    return r1.chr == r2.chr && r1.start == r2.start && r1.end == r2.end;
}

sv_t parse_sv(string& line) {
    char id[100];
    char chr1[100], chr2[100];
    int pos1, pos2;
    char c;
    sscanf(line.data(), "%s %s %d %c %s %d %c", id, chr1, &pos1, &c, chr2, &pos2, &c);
    return sv_t(id, bp_t(chr1,pos1), bp_t(chr2,pos2));
}

repeat_t parse_rep(string& line) {
    char chr[100];
    int start, end;
    sscanf(line.data(), "%s %s %d %d", chr, chr, &start, &end);
    return repeat_t(chr, start, end);
}

bool in_rep(sv_t sv, repeat_t r) {
    return sv.bp1.chr == r.chr && r.start <= sv.bp1.pos && r.end >= sv.bp2.pos;
}

bool intersect_rep(sv_t sv, repeat_t r) {
    return sv.bp1.chr == r.chr && sv.bp1.pos <= r.end && sv.bp2.pos >= r.start;
}

string print_sv(sv_t& sv) {
    stringstream ss;
    ss << sv.bp1.chr << ":" << sv.bp1.pos << "-" << sv.bp2.pos;
    return ss.str();
}
string print_rep(repeat_t& r) {
    stringstream ss;
    ss << r.chr << ":" << r.start << "-" << r.end;
    return ss.str();
}

int dist(sv_t& sv1, sv_t& sv2) {
    if (sv1.bp1.chr != sv2.bp1.chr) return INT_MAX;
    return abs(sv1.bp1.pos-sv2.bp1.pos) + abs(sv1.bp2.pos-sv2.bp2.pos);
}

int get_overlap(sv_t& sv1, sv_t& sv2) {
    return max(0, min(sv1.bp2.pos, sv2.bp2.pos)-max(sv1.bp1.pos, sv2.bp1.pos));
}

int main(int argc, char* argv[]) {

    if (argc < 6) {
        cout << "Given a SV file with benchmark deletions and one with the called ones, reports for each benchmark deletions if it is being called." << endl;
        cout << "Usage: exec benchmark_file called_file ucsc_repeats ref_genome maxdist" << endl;
        return 0;
    }

    ifstream benchmark_f(argv[1]), called_f(argv[2]), rep_f(argv[3]);
    int maxdist = stoi(argv[5]);
    double overlap_frac = 0.9;
    if (argc == 7) overlap_frac = std::stod(argv[6]);

    unordered_map<string, string> chrs;
    
    FILE* fastaf = fopen(argv[4], "r");
    kseq_t *seq = kseq_init(fileno(fastaf));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        chrs[string(seq->name.s)] = string(seq->seq.s);
    }
    kseq_destroy(seq);

    vector<sv_t> benchmark_svs;
    unordered_map<string, vector<sv_t> > called_svs;

    string line;
    while (getline(benchmark_f, line)) {
        if (line.find(" DEL") == string::npos) continue;
        benchmark_svs.push_back(parse_sv(line));
    }
    while (getline(called_f, line)) {
        if (line.find(" DEL") == string::npos) continue;
        sv_t del = parse_sv(line);
        called_svs[del.bp1.chr].push_back(del);
    }

    unordered_map<string, vector<repeat_t>> reps;
    unordered_map<string, vector<Interval<repeat_t>>> reps_iv;
    unordered_map<string, IntervalTree<repeat_t>*> reps_i;
    while (getline(rep_f, line)) {
        if (line[0] == '#') continue;
        repeat_t r = parse_rep(line);
        reps[r.chr].push_back(r);
        reps_iv[r.chr].push_back(Interval<repeat_t>(r.start, r.end, r));
    }
    
    for (auto it = reps.begin(); it != reps.end(); it++) {
        reps_i[it->first] = new IntervalTree<repeat_t>(reps_iv[it->first]);
    }
    
    Aligner aligner(2,2,4,1);
    Filter filter(false, false, 0, 32767);
    Alignment alignment;

    for (int i = 0; i < benchmark_svs.size(); i++) {
        sv_t bsv = benchmark_svs[i];
        
        bool skip = false;
        for (int j = 0; j < called_svs[bsv.bp1.chr].size(); j++) {
            sv_t csv = called_svs[bsv.bp1.chr][j];
            if (dist(bsv,csv) <= maxdist && get_overlap(bsv,csv) >= min(bsv.size(), csv.size())*overlap_frac) {
                cout << bsv.id << " " << csv.id << endl;
                skip = true;
               // break;
            }
        }
        
        if (reps_i[bsv.bp1.chr] == NULL) {
            if (!skip) cout << bsv.id << " NONE" << endl;
            continue;
        }

        vector<Interval<repeat_t>> intervals_temp = reps_i[bsv.bp1.chr]->findOverlapping(bsv.bp1.pos, bsv.bp2.pos);
        vector<repeat_t> reps_containing;
        for (int j = 0; j < intervals_temp.size(); j++) {
            repeat_t r = intervals_temp[j].value;
            if (in_rep(bsv, r)) {
                reps_containing.push_back(r);
            }
        }
        
        for (int j = 0; j < called_svs[bsv.bp1.chr].size(); j++) {
            sv_t csv = called_svs[bsv.bp1.chr][j];
            if (abs(bsv.size()-csv.size()) > maxdist) continue;

            for (int k = 0; k < reps_containing.size(); k++) {
                repeat_t r = reps_containing[k];
                if (intersect_rep(csv, r)) {
                    string bseq = chrs[bsv.bp1.chr].substr(bsv.bp1.pos, bsv.bp2.pos-bsv.bp1.pos);
                    string cseq = chrs[csv.bp1.chr].substr(csv.bp1.pos, csv.bp2.pos-csv.bp1.pos);
                    aligner.Align(bseq.data(), cseq.data(), cseq.length(), filter, &alignment);
                    if (alignment.sw_score >= min(bseq.length(), cseq.length())) {
                        //cout << print_sv(bsv) << " " << print_sv(csv) << " REP ";
                        //cout << min(bseq.length(), cseq.length()) << " " << alignment.sw_score << endl;
                        cout << bsv.id << " " << csv.id << " REP" << endl;
                        skip = true;
                        break;
                    }
                }
            }

            //if (skip) break;
        }
        
        if (skip) continue;
        
        cout << bsv.id << " NONE" << endl;
    }
}
