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

const int INF = 1000000;

unordered_map<string, string> chrs;
int maxdist;
Aligner aligner(2,2,4,1);
Filter filter(false, false, 0, 32767);
Alignment alignment;

struct bp_t {
	string chr;
	int pos;

	bp_t(string chr, int pos) : chr(chr), pos(pos) {}
};

struct sv_t {
	string id;
	bp_t bp1, bp2;
    string svtype, seq;

	sv_t(string id, bp_t bp1, bp_t bp2, string svtype, string seq) : id(id), bp1(bp1), bp2(bp2), svtype(svtype), seq(seq) {}

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
    char svtype[100], seq[1000000];
	sscanf(line.data(), "%s %s %d %c %s %d %c %s %s", id, chr1, &pos1, &c, chr2, &pos2, &c, &svtype, &seq);
	return sv_t(id, bp_t(chr1,pos1), bp_t(chr2,pos2), svtype, seq);
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

int dist(bp_t& bp1, bp_t& bp2) {
    if (bp1.chr != bp2.chr) return INF;
    return abs(bp1.pos - bp2.pos);
}
int dist(sv_t& sv1, sv_t& sv2) {
	if (sv1.bp1.chr != sv2.bp1.chr) return INF;
	return dist(sv1.bp1, sv2.bp1) + dist(sv1.bp2, sv2.bp2);
}

string reverse_comple(string& seq) {
    char rc[1000000];
    int32_t end = seq.length(), start = 0;
    static const int8_t rc_table[128] = {
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4,  4, 4,  4,  4,  4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
        4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 84, 4, 71, 4,  4,  4, 67, 4, 4, 4, 4,  4, 4, 4, 4,
        4, 4,  4, 4,  65, 65, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
    };
    rc[end] = '\0';
    -- end;
    while (start < end) {
        rc[start] = (char)rc_table[(int8_t)seq[end]];
        rc[end] = (char)rc_table[(int8_t)seq[start]];
        ++ start;
        -- end;
    }
    if (start == end) rc[start] = (char)rc_table[(int8_t)seq[start]];
    return string(rc);
}

bool check_regions(sv_t& bsv, sv_t& csv, bool first_is_stable = false) {
    if (csv.svtype == "INS" || csv.svtype == "NOV" || bsv.size() > 50000 || csv.size() > 50000) return true;
    
    string bseq;
    if (bsv.svtype == "DUP") {
        bseq = chrs[bsv.bp1.chr].substr(bsv.bp1.pos, bsv.bp2.pos-bsv.bp1.pos);
    }

    string cseq;
    if (csv.bp1.pos >= csv.bp2.pos) return false;
    if (csv.svtype == "DUP") {
        cseq = chrs[csv.bp1.chr].substr(csv.bp1.pos, csv.bp2.pos-csv.bp1.pos);
    } else if (csv.svtype == "TRA") {
        // the stable bp in a TRA is the one near the insertion point, unstable the other one
        bp_t& unstable_bp = first_is_stable ? csv.bp2 : csv.bp1;
        int a = unstable_bp.pos - maxdist - bsv.seq.length();
        a = max(a, 0);
        int b = unstable_bp.pos + maxdist + bsv.seq.length();
        cseq = chrs[unstable_bp.chr].substr(a, b-a);
    }
    
    bseq = bseq + bseq;
    aligner.Align(bseq.data(), cseq.data(), cseq.length(), filter, &alignment);
    if (alignment.sw_score >= min(bsv.seq.length(), cseq.length())) return true;

    cseq = reverse_comple(cseq);
    aligner.Align(bseq.data(), cseq.data(), cseq.length(), filter, &alignment);
    return alignment.sw_score >= min(bsv.seq.length(), cseq.length());
}

int main(int argc, char* argv[]) {
    
	if (argc < 6) {
		cout << "Given a SV file with benchmark insertions and one with the called ones, reports for each benchmark insertion if it is being called." << endl;
		cout << "Usage: exec benchmark_file called_file ucsc_repeats ref_genome maxdist" << endl;
		return 0;
	}

	ifstream benchmark_f(argv[1]), called_f(argv[2]), rep_f(argv[3]);
	maxdist = stoi(argv[5]);

    FILE* fastaf = fopen(argv[4], "r");
    kseq_t *seq = kseq_init(fileno(fastaf));
    int l;
    while ((l = kseq_read(seq)) >= 0) {
        chrs[string(seq->name.s)] = string(seq->seq.s);
    }
    kseq_destroy(seq);

	vector<sv_t> benchmark_svs, called_svs;

	string line;
	while (getline(benchmark_f, line)) {
		benchmark_svs.push_back(parse_sv(line));
	}
	while (getline(called_f, line)) {
		if (line.find(" INS") == string::npos && line.find(" DUP") == string::npos && line.find(" TRA") == string::npos) continue;
		called_svs.push_back(parse_sv(line));
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

	for (int i = 0; i < benchmark_svs.size(); i++) {
		sv_t bsv = benchmark_svs[i];

        bool skip = false;
        for (int j = 0; j < called_svs.size(); j++) {
            sv_t csv = called_svs[j];
            if (((csv.svtype == "INS") && dist(bsv,csv) <= maxdist) ||
                ((csv.svtype == "DUP" || csv.svtype == "TRA") &&
                min(dist(csv.bp1, bsv.bp1), dist(csv.bp2, bsv.bp2)) <= maxdist)) {
                if (check_regions(bsv, csv, dist(csv.bp1, bsv.bp1) <= dist(csv.bp2, bsv.bp2))) {
                    cout << bsv.id << " " << csv.id << endl;
                    skip = true;
                }
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

        for (int j = 0; j < called_svs.size(); j++) {
            sv_t csv = called_svs[j];
            if (bsv.size() < csv.size()-maxdist || bsv.bp1.chr != csv.bp1.chr) continue;

            for (int k = 0; k < reps_containing.size(); k++) {
                repeat_t r = reps_containing[k];
                if (intersect_rep(csv, r)) {
                    if (check_regions(bsv, csv, dist(csv.bp1, bsv.bp1) <= dist(csv.bp2, bsv.bp2))) {
                        cout << bsv.id << " " << csv.id << " REP" << endl;
                        skip = true;
                        break;
                    }
				}
			}

//            if (skip) break;
		}
        
        if (skip) continue;
        
        cout << bsv.id << " NONE" << endl;
	}
}
