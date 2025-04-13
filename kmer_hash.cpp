#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <vector>
#include <fstream>
#include <upcxx/upcxx.hpp>
#include <upcxx/dist_object.hpp>
#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

int main(int argc, char** argv) {
    upcxx::init();
    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }
    std::string kmer_fname = argv[1];
    std::string run_type = "";
    if (argc >= 3) {
        run_type = argv[2];
    }
    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = argv[3];
    }

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " has " + std::to_string(ks) +
                                 "-mers, but this binary expects " + std::to_string(KMER_LEN) + "-mers.");
    }

    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = static_cast<size_t>(n_kmers * (1.0 / 0.5));
    HashMap hashmap(hash_table_size);

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());

    upcxx::barrier();
    auto start_insert = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> start_nodes;
    for (auto &kmer : kmers) {
        bool ok = hashmap.insert(kmer);
        if (!ok) {
            throw std::runtime_error("Error: HashMap is full!");
        }
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }

    upcxx::barrier();
    auto end_insert = std::chrono::high_resolution_clock::now();
    double insert_time = std::chrono::duration<double>(end_insert - start_insert).count();
    if (run_type != "test" && upcxx::rank_me() == 0) {
        BUtil::print("Finished inserting in %lf\n", insert_time);
    }

    upcxx::dist_object<std::vector<kmer_pair>> dstart(start_nodes);
    upcxx::barrier();
    auto start_read = std::chrono::high_resolution_clock::now();

    std::vector<kmer_pair> merged_starts;
    if (upcxx::rank_me() == 0) {
        merged_starts = dstart.fetch(0).wait();
        for (int r = 1; r < upcxx::rank_n(); r++) {
            auto remote_vec = dstart.fetch(r).wait();
            merged_starts.insert(merged_starts.end(), remote_vec.begin(), remote_vec.end());
        }
    }

    std::list<std::list<kmer_pair>> contigs;
    if (upcxx::rank_me() == 0) {
        for (auto &start_k : merged_starts) {
            std::list<kmer_pair> contig;
            contig.push_back(start_k);
            while (contig.back().forwardExt() != 'F') {
                kmer_pair next_k;
                bool found = hashmap.find(contig.back().next_kmer(), next_k);
                if (!found) {
                    throw std::runtime_error("Error: k-mer not found in hashmap.");
                }
                contig.push_back(next_k);
            }
            contigs.push_back(contig);
        }
    }

    upcxx::barrier();
    auto end_read = std::chrono::high_resolution_clock::now();
    auto total_end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> read_time = end_read - start_read;
    std::chrono::duration<double> total_time = total_end - start_insert;

    if (upcxx::rank_me() == 0) {
        int numKmers = 0;
        for (auto &c : contigs) {
            numKmers += (int)c.size();
        }
        if (run_type != "test") {
            BUtil::print("Assembled in %lf total\n", total_time.count());
        }
        if (run_type == "verbose") {
            printf("Rank 0: reconstructed %lu contigs with %d nodes. (read %lf, insert %lf, total %lf)\n",
                   contigs.size(), numKmers, read_time.count(), insert_time, total_time.count());
        }
        if (run_type == "test") {
            std::ofstream fout(test_prefix + "_0.dat");
            for (auto &c : contigs) {
                fout << extract_contig(c) << "\n";
            }
            fout.close();
        }
    }

    upcxx::finalize();
    return 0;
}
