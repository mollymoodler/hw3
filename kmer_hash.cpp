#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <list>
#include <numeric>
#include <set>
#include <upcxx/upcxx.hpp>
#include <vector>
#include <fstream>
#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"
#include "butil.hpp"

int main(int argc, char **argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = argv[1];
    std::string run_type = (argc >= 3 ? argv[2] : "");
    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = argv[3];
    }

    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) + "-mers. Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = n_kmers * 2; // 50% load factor -> table size = 2 * n_kmers
    HashMap hashmap(hash_table_size);

    if (run_type == "verbose") {
        BUtil::print("Initializing hash table of size %zu for %zu kmers.\n", hash_table_size, n_kmers);
    }

    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    if (run_type == "verbose") {
        BUtil::print("Finished reading kmers.\n");
    }

    auto start = std::chrono::high_resolution_clock::now();
    std::vector<kmer_pair> start_nodes;
    for (auto &kmer : kmers) {
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap is full!");
        }
        if (kmer.backwardExt() == 'F') {
            start_nodes.push_back(kmer);
        }
    }
    auto end_insert = std::chrono::high_resolution_clock::now();
    upcxx::barrier();

    double insert_time = std::chrono::duration<double>(end_insert - start).count();
    if (run_type != "test")
        BUtil::print("Finished inserting in %lf seconds\n", insert_time);
    upcxx::barrier();

    auto start_read = std::chrono::high_resolution_clock::now();
    std::list<std::list<kmer_pair>> contigs;
    for (const auto &start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);
        while (contig.back().forwardExt() != 'F') {
            kmer_pair kmer;
            bool found = hashmap.find(contig.back().next_kmer(), kmer);
            if (!found) {
                throw std::runtime_error("Error: k-mer not found in hashmap.");
            }
            contig.push_back(kmer);
        }
        contigs.push_back(contig);
    }
    auto end_read = std::chrono::high_resolution_clock::now();
    upcxx::barrier();
    auto total = std::chrono::high_resolution_clock::now() - start;
    double read_time = std::chrono::duration<double>(end_read - start_read).count();

    int numKmers = std::accumulate(
        contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair> &contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Assembled in %lf seconds total.\n",
                     std::chrono::duration<double>(total).count());
    }

    if (run_type == "verbose") {
        printf("Rank %d reconstructed %zu contigs with %d nodes from %zu start nodes. (read: %lf, insert: %lf, total: %lf)\n",
               upcxx::rank_me(), contigs.size(), numKmers, start_nodes.size(),
               read_time, std::chrono::duration<double>(end_insert - start).count(),
               std::chrono::duration<double>(total).count());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto &contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}