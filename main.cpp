#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <sstream>
#include <vector>
#include "SubstitutionMatrix.h"
#include "FwBwAligner.h"

using namespace std;


unordered_map<string, string> load_sequences(const string& filename) {
    unordered_map<string, string> sequences;
    ifstream file(filename);
    string sid = "";
    string line;
    while (getline(file, line)) {
        line.erase(0, line.find_first_not_of(">"));
        line.erase(line.find_last_not_of('\n') + 1);
        if (sid.empty()) {
            sid = line;
        } else {
            sequences[sid] = line;
            sid = "";
        }
    }
    return sequences;
}

std::string format_alignment(
    const FwBwAligner::s_align& result,
    const std::string& query,
    const std::string& target,
    const std::string& queryName,
    const std::string& targetName
) {
    std::string aln_q = "";
    std::string aln_t = "";
    aln_q.reserve(result.cigarLen);
    aln_t.reserve(result.cigarLen);
    size_t idx_i = result.qStartPos1;
    size_t idx_j = result.dbStartPos1;
    for (const char& c : result.cigar) {
        if (c == 'M') {
            aln_q.push_back(query[idx_i++]);
            aln_t.push_back(target[idx_j++]);
        } else if (c == 'D') {
            aln_q.push_back(query[idx_i++]);
            aln_t.push_back('-');
        } else {
            aln_q.push_back('-');
            aln_t.push_back(target[idx_j++]);
        }
    }
    std::string fasta = ">";
    fasta.append(queryName);
    fasta.append(1, '\n');
    fasta.append(aln_q);
    fasta.append("\n>");
    fasta.append(targetName);
    fasta.append(1, '\n');
    fasta.append(aln_t);
    fasta.append(1, '\n');

    return fasta;
}


int main() {
    unordered_map<string, string> sid2seq_ss = load_sequences("./data/ss.fasta");
    unordered_map<string, string> sid2seq_aa = load_sequences("./data/aa.fasta");

    string query = "d1ufaa2"; //d1ufaa2
    string target = "d2c1ia1"; //d2c1ia1
    // string query = "d1ufaa2";
    // string target = "d2c1ia1";

    string aa1 = sid2seq_aa[query];
    string aa2 = sid2seq_aa[target];
    string ss1 = sid2seq_ss[query];
    string ss2 = sid2seq_ss[target];

    size_t rows = ss1.size();
    size_t cols = ss2.size();
    std::cout << "rows: " << rows << " cols: " << cols << std::endl;

    SubstitutionMatrix subMat3Di("./data/mat3di.out", 2.1, 0.0);
    SubstitutionMatrix subMat("./data/blosum62.out", 1.4, 0.0);

    FwBwAligner aligner(rows, cols, subMat3Di, subMat);
    FwBwAligner::s_align result = aligner.align(ss1, aa1, ss2, aa2);

    std::string fasta = format_alignment(result, aa1, aa2, query, target);
    std::cout << fasta;

    return 0;
}
