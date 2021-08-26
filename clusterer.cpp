#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <set>

struct sv_t {
    std::string id, chr, type, seq, sample;
    int start, end;
    char str1, str2;

    sv_t(std::string id, std::string chr, int start, char str1, int end, char str2, std::string type, std::string seq,
         std::string sample)
    : id(id), chr(chr), start(start), str1(str1), end(end), str2(str2), type(type), seq(seq), sample(sample) {}

    std::string to_string() {
        std::stringstream ss;
        ss << sample << "." << id << "," << chr << ":" << start << "-" << end;
        return ss.str();
    }
};
bool operator < (const sv_t& sv1, const sv_t& sv2) {
    if (sv1.chr != sv2.chr) return sv1.chr < sv2.chr;
    return sv1.start < sv2.start;
}
int distance(const sv_t& sv1, const sv_t& sv2) {
    if (sv1.chr != sv2.chr) return INT32_MAX;
    return abs(sv1.start-sv2.start) + abs(sv1.end-sv2.end);
}

struct idx_size_t {
    int idx, size;

    idx_size_t() : idx(-1), size(0) {}
    idx_size_t(int idx, int size) : idx(idx), size(size) {}
};
bool operator > (const idx_size_t& s1, const idx_size_t& s2) {
    return s1.size > s2.size;
};

std::vector<sv_t> svs;
std::vector<std::vector<int>> compatibility_list;
int max_dist, min_overlap;

int overlap(const sv_t& sv1, const sv_t& sv2) {
    return std::max(0, std::min(sv1.end, sv2.end)-std::max(sv1.start, sv2.start));
}

bool is_compatible(const sv_t& sv1, const sv_t& sv2) {
    return distance(sv1, sv2) <= max_dist && overlap(sv1, sv2) >= min_overlap && abs((int)(sv1.seq.length()-sv2.seq.length())) < 100;
}

int median(std::vector<int>& v) {
    if (v.size()%2 == 0) {
        return (v[v.size()/2-1] + v[v.size()/2])/2;
    } else {
        return v[v.size()/2];
    }
}

std::vector<std::vector<int>> compute_minimal_clique_cover(int start, int end) {
    std::vector<idx_size_t> idx_by_neighborhood_size;
    for (int i = start; i < end; i++) {
        idx_by_neighborhood_size.emplace_back(i, compatibility_list[i].size());
    }
    std::sort(idx_by_neighborhood_size.begin(), idx_by_neighborhood_size.end(), std::greater<idx_size_t>());

    std::vector<std::vector<int>> cliques;
    int cliques_counter = 0;
    std::vector<int> cliques_idx(svs.size(), -1);
    for (idx_size_t& curr_idx : idx_by_neighborhood_size) {
        // find cliqued neighbors
        std::vector<idx_size_t> cliqued_neighbors;
        for (int adj_idx : compatibility_list[curr_idx.idx]) {
            int clique_idx = cliques_idx[adj_idx];
            if (clique_idx != -1 && is_compatible(svs[curr_idx.idx], svs[adj_idx])) {
                cliqued_neighbors.emplace_back(adj_idx, cliques[clique_idx].size());
            }
        }

        // check if I can join any neighboring cliques
        std::sort(cliqued_neighbors.begin(), cliqued_neighbors.end(), std::greater<idx_size_t>());
        bool used = false;
        for (idx_size_t& neighbor_idx : cliqued_neighbors) {
            int neighboor_clique_idx = cliques_idx[neighbor_idx.idx];
            std::vector<int> &neighboor_clique = cliques[neighboor_clique_idx];
            bool belongs_to_clique = true;
            for (int clique_elem_idx : neighboor_clique) {
                if (!is_compatible(svs[curr_idx.idx], svs[clique_elem_idx])) {
                    belongs_to_clique = false;
                    break;
                }
            }

            if (belongs_to_clique) {
                neighboor_clique.push_back(curr_idx.idx);
                used = true;
                break;
            }
        }

        if (!used) {
            cliques.push_back(std::vector<int>());
            cliques[cliques_counter].push_back(curr_idx.idx);
            cliques_idx[curr_idx.idx] = cliques_counter;
            cliques_counter++;
        }
    }

    return cliques;
}

int cluster_id = 0;
void print_cliques(std::vector<std::vector<int>>& cliques) {
    for (auto& clique : cliques) {
        std::vector<int> starts, ends, seq_len;
        std::set<std::string> unique_samples;
        for (int idx : clique) {
            starts.push_back(svs[idx].start);
            ends.push_back(svs[idx].end);
            seq_len.push_back(svs[idx].seq.length());
            unique_samples.insert(svs[idx].sample);
        }
        std::sort(starts.begin(), starts.end());
        std::sort(ends.begin(), ends.end());
        std::cout << "CLUSTER_" << cluster_id++ << " " << svs[clique[0]].chr << " " << median(starts) << " " << svs[clique[0]].str1;
        std::cout << " " << svs[clique[0]].chr << " " << median(ends) << " " << svs[clique[0]].str2 << " ";
        std::cout << svs[clique[0]].type << " " << unique_samples.size() << " " << clique.size() << " " << median(seq_len) << " ";
        std::vector<sv_t> clique_svs;
        for (int idx : clique) {
            clique_svs.push_back(svs[idx]);
        }
        std::sort(clique_svs.begin(), clique_svs.end());
        for (sv_t sv : clique_svs) {
            std::cout << " " << sv.to_string();
        }
        std::cout << std::endl;
    }
}

int main(int argc, char* argv[]) {
    std::ifstream sv_fin(argv[1]);
    max_dist = std::stoi(argv[2]);
    min_overlap = std::stoi(argv[3]);

    std::string line;
    char id[128], chr[1024], type[16], sample[128], precise[100], seq[100000];
    char str1, str2;
    int start, end;
    while (getline(sv_fin, line)) {
        sscanf(line.c_str(), "%s %s %d %c %*s %d %c %s %s %s %s", id, chr, &start, &str1, &end, &str2, type, seq, /*precise,*/
               sample);
        svs.emplace_back(id, chr, start, str1, end, str2, type, seq, sample);
        compatibility_list.push_back(std::vector<int>());
    }

    std::sort(svs.begin(), svs.end());
    int i_prev = 0;
    for (int i = 0; i < svs.size(); i++) {
        sv_t& sv1 = svs[i];
        for (int j = i+1; j < svs.size(); j++) {
            sv_t& sv2 = svs[j];
            if (sv1.chr != sv2.chr || sv2.start-sv1.start > max_dist) break;
            if (is_compatible(sv1, sv2)) {
                compatibility_list[i].push_back(j);
                compatibility_list[j].push_back(i);
            }
        }
        compatibility_list[i].shrink_to_fit();

        if (i > 0 && (sv1.chr != svs[i-1].chr || sv1.start-svs[i-1].start > max_dist)) {
            std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, i);
            print_cliques(cliques);

            for (int j = i_prev; j < i; j++) {
                compatibility_list[j].clear();
                compatibility_list[j].shrink_to_fit();
            }
            i_prev = i;
        }
    }
    std::vector<std::vector<int>> cliques = compute_minimal_clique_cover(i_prev, svs.size());
    print_cliques(cliques);
}
