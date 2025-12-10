#include <bits/stdc++.h>
using namespace std;

// bedavg.cpp
// Fast, low-memory C++ program to compute per-region averages for many samples.
// Usage:
//   bedavg met.tsv regions.bed out.tsv [threads] [--bedzero|--1based] [--inclusive]
// - met.tsv: tab-separated with header: chr\tpos\tSample1\tSample2\t...
// - regions.bed: BED-like: chr\tstart\tend\tname (start,end integers)
// Defaults assume met is sorted by chr then pos. regions.bed may be in any order.
// By default: BED is treated as 0-based half-open (start <= pos-1 < end).
// If --1based is provided, intervals are treated as 1-based inclusive (start <= pos <= end).
// If --inclusive is provided, the check will be start <= pos <= end regardless of --bedzero/--1based.

struct Interval {
    long long start;
    long long end;
    string name;
    int index; // position in vector
};

// Helper: split by tab (fast)
static inline void split_tab(const string &s, vector<string> &out) {
    out.clear();
    size_t i = 0, n = s.size();
    while (i < n) {
        size_t j = i;
        while (j < n && s[j] != '\t' && s[j] != '\r' && s[j] != '\n') ++j;
        out.emplace_back(s.data() + i, j - i);
        // skip separators
        i = j + 1;
        while (i < n && (s[i] == '\r' || s[i] == '\n')) ++i;
    }
}

// Trim helper
static inline string trim(const string &s) {
    size_t a = 0, b = s.size();
    while (a < b && isspace((unsigned char)s[a])) ++a;
    while (b > a && isspace((unsigned char)s[b-1])) --b;
    return s.substr(a, b-a);
}

// Read bed file and group intervals by chromosome
void read_bed(const string &bedfile, unordered_map<string, vector<Interval>> &bed_map) {
    ifstream in(bedfile);
    if (!in) { cerr << "Error opening bed file: " << bedfile << "\n"; exit(1); }
    string line;
    vector<string> f;
    long long lineno = 0;
    while (getline(in, line)) {
        ++lineno;
        if (line.empty() || line[0] == '#') continue;
        split_tab(line, f);
        if (f.size() < 3) continue;
        string chr = f[0];
        long long start = stoll(trim(f[1]));
        long long end = stoll(trim(f[2]));
        // string name = (f.size() >= 4 ? f[3] : (chr + ":" + to_string(start) + "-" + to_string(end)));
        string name = chr + "-" + to_string(start) + "-" + to_string(end);
        Interval iv{start, end, name, 0};
        bed_map[chr].push_back(iv);
    }
    // assign indices
    for (auto &p : bed_map) {
        auto &v = p.second;
        sort(v.begin(), v.end(), [](const Interval &a, const Interval &b){
            if (a.start != b.start) return a.start < b.start;
            return a.end < b.end;
        });
        for (size_t i=0;i<v.size();++i) v[i].index = (int)i;
    }
}

// Determine byte offsets in met file for each chromosome block (assumes met is sorted by chr)
void index_met_by_chr(const string &metfile, unordered_map<string, pair<long long,long long>> &chr_offsets, string &header) {
    ifstream in(metfile, ios::binary);
    if (!in) { cerr << "Error opening met file: " << metfile << "\n"; exit(1); }
    // read header line
    long long pos0 = in.tellg();
    if (!getline(in, header)) { cerr << "Empty met file or can't read header\n"; exit(1); }
    long long after_header = in.tellg();
    string line;
    vector<string> f;
    string current_chr = "";
    long long block_start = -1;
    while (getline(in, line)) {
        if (line.empty()) continue;
        split_tab(line, f);
        if (f.size() < 2) continue;
        string chr = f[0];
        if (chr != current_chr) {
            if (!current_chr.empty()) {
                chr_offsets[current_chr].second = in.tellg(); // position after reading line; approximate
            }
            current_chr = chr;
            chr_offsets[chr].first = in.tellg() - line.size() - 1; // adjust to include first row
            // The recorded .first is actually after the first line of that chr; okay for seeking
        }
    }
    // last chromosome end is EOF (we'll represent with -1 meaning till EOF)
    if (!current_chr.empty()) chr_offsets[current_chr].second = -1;
}

// Worker for one chromosome
void process_chromosome(const string &metfile, const string &chr, const vector<Interval> &intervals_in,
                        long long met_start_offset, long long met_end_offset, const vector<string> &sample_names,
                        bool bed_zero_based, bool inclusive, ostream &out, mutex &out_mtx) {
    if (intervals_in.empty()) return;
    // Copy intervals (we will store sums/counts)
    size_t nreg = intervals_in.size();
    struct RegAccum { vector<double> sums; long long count; string name; long long start,end; };
    vector<RegAccum> regs;
    regs.reserve(nreg);
    for (const auto &iv : intervals_in) {
        RegAccum r; r.sums.assign(sample_names.size(), 0.0); r.count = 0; r.name = iv.name; r.start = iv.start; r.end = iv.end; regs.push_back(move(r));
    }
    // active intervals: store indices
    vector<int> active;
    size_t next_interval = 0; // index in intervals_in to add when start <= pos

    ifstream in(metfile);
    if (!in) { cerr << "Error opening met file in thread for chr "<<chr<<"\n"; return; }
    if (met_start_offset > 0) in.seekg(met_start_offset);
    string line;
    vector<string> f;

    auto pos_in_interval = [&](long long pos, long long istart, long long iend)->bool{
        if (inclusive) return (istart <= pos && pos <= iend);
        if (bed_zero_based) {
            // bed 0-based half-open -> match if istart <= pos-1 < iend --> pos-1 >= istart && pos-1 < iend
            long long p0 = pos - 1;
            return (istart <= p0 && p0 < iend);
        } else {
            // assume 1-based inclusive start,end
            return (istart <= pos && pos <= iend);
        }
    };

    // read lines until met_end_offset (if -1 then till EOF)
    while (true) {
        long long here = in.tellg();
        if (met_end_offset != -1 && here >= met_end_offset) break;
        if (!getline(in, line)) break;
        if (line.empty()) continue;
        split_tab(line, f);
        if (f.size() < 2) continue;
        string rchr = f[0];
        if (rchr != chr) break; // reached next chr
        long long pos = stoll(f[1]);
        // advance next_interval to include intervals whose start <= pos (accounting for bed_zero/1based)
        while (next_interval < intervals_in.size()) {
            long long ivstart = intervals_in[next_interval].start;
            long long testpos = pos;
            // We need to check if ivstart <= pos (using a conservative criterion):
            // For bed_zero_based and non-inclusive default, we add interval when ivstart <= pos-1 => ivstart <= pos-1
            bool should_add = false;
            if (inclusive) should_add = (ivstart <= pos);
            else if (bed_zero_based) should_add = (ivstart <= pos-1);
            else should_add = (ivstart <= pos);
            if (should_add) {
                active.push_back((int)next_interval);
                ++next_interval;
            } else break;
        }
        // remove expired from active (those where pos no longer in interval)
        for (size_t i=0;i<active.size();) {
            int idx = active[i];
            long long iend = intervals_in[idx].end;
            bool still = pos_in_interval(pos, intervals_in[idx].start, iend);
            if (!still) {
                // remove by swap
                active[i] = active.back(); active.pop_back();
            } else ++i;
        }
        if (active.empty()) continue;
        // parse sample values starting at column 2
        size_t nsamples = sample_names.size();
        // Defensive: if columns fewer than expected, skip
        if (f.size() < 2 + nsamples) continue;
        for (int ai : active) {
            auto &acc = regs[ai];
            for (size_t si = 0; si < nsamples; ++si) {
                const string &sv = f[2+si];
                if (sv.empty() || sv == "." || sv == "-" || sv == "NA" || sv == "nan") continue;
                double v = atof(sv.c_str());
                acc.sums[si] += v;
            }
            acc.count += 1;
        }
    }

    // write results to out with mutex
    {
        lock_guard<mutex> lg(out_mtx);
        // header for first writer? We'll assume main program writes header. Here write rows
        for (size_t i=0;i<regs.size();++i) {
            auto &r = regs[i];
            out << r.name;
            for (size_t si=0; si<sample_names.size(); ++si) {
                out << '\t';
                if (r.count == 0) out << "NA";
                else {
                    double avg = r.sums[si] / (double)r.count;
                    // print with 6 decimal places
                    out.setf(std::ios::fixed); out.precision(6);
                    out << avg;
                    out.unsetf(std::ios::floatfield);
                }
            }
            out << '\n';
        }
    }
}

int main(int argc, char** argv) {
    if (argc < 4) {
        cerr << "Usage: "<<argv[0]<<" met.tsv regions.bed out.tsv [threads] [--bedzero|--1based] [--inclusive]\n";
        return 1;
    }
    string metfile = argv[1];
    string bedfile = argv[2];
    string outfile = argv[3];
    int nthreads = 1;
    bool bed_zero_based = true;
    bool inclusive = false;
    for (int i=4;i<argc;++i) {
        string a = argv[i];
        if (a == "--1based") { bed_zero_based = false; }
        else if (a == "--bedzero") { bed_zero_based = true; }
        else if (a == "--inclusive") { inclusive = true; }
        else {
            // if numeric, set threads
            bool isnumeric = true; for (char c: a) if (!isdigit((unsigned char)c)) isnumeric=false;
            if (isnumeric) nthreads = stoi(a);
        }
    }

    // Read bed
    unordered_map<string, vector<Interval>> bed_map;
    read_bed(bedfile, bed_map);
    if (bed_map.empty()) { cerr << "No intervals read from "<<bedfile<<"\n"; return 1; }

    // Index met file by chromosome offsets
    unordered_map<string, pair<long long,long long>> chr_offsets; // chr -> (start_offset, end_offset or -1)
    string header;
    index_met_by_chr(metfile, chr_offsets, header);
    if (header.empty()) { cerr << "Empty header in met file"<<"\n"; return 1; }
    // parse header for sample names
    vector<string> head_fields; split_tab(header, head_fields);
    if (head_fields.size() < 3) { cerr << "met header must have at least chr,pos,sample...\n"; return 1; }
    vector<string> sample_names;
    for (size_t i=2;i<head_fields.size();++i) sample_names.push_back(head_fields[i]);

    // Prepare chromosomes to process: intersection of bed_map keys and chr_offsets keys
    vector<string> chroms;
    for (auto &p : bed_map) {
        if (chr_offsets.find(p.first) != chr_offsets.end()) chroms.push_back(p.first);
    }
    if (chroms.empty()) { cerr << "No chromosomes in common between met and bed\n"; return 1; }

    // open output and write header
    ofstream out(outfile);
    if (!out) { cerr << "Cannot open output file "<<outfile<<"\n"; return 1; }
    out << "name";
    for (auto &s: sample_names) out << '\t' << s;
    out << '\n';

    // Thread pool
    mutex out_mtx;
    vector<thread> pool;
    atomic<size_t> next_chr_idx(0);

    auto worker = [&]() {
        while (true) {
            size_t idx = next_chr_idx.fetch_add(1);
            if (idx >= chroms.size()) break;
            string chr = chroms[idx];
            auto &intervals = bed_map[chr];
            long long moff_s = chr_offsets[chr].first;
            long long moff_e = chr_offsets[chr].second;
            process_chromosome(metfile, chr, intervals, moff_s, moff_e, sample_names, bed_zero_based, inclusive, out, out_mtx);
        }
    };

    int use_threads = max(1, min((int)chroms.size(), nthreads));
    for (int t=0;t<use_threads;++t) pool.emplace_back(worker);
    for (auto &th: pool) th.join();

    out.close();
    cerr << "Done. Results written to "<<outfile<<"\n";
    return 0;
}
