#include <chrono>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <list>
#include <numeric>
#include <set>
#include <vector>
#include <stdexcept>
#include <upcxx/upcxx.hpp>
#include "hash_map.hpp"
#include "kmer_t.hpp"
#include "read_kmers.hpp"  
#include "butil.hpp"

#define MAX_CHAIN_LENGTH 1000  // Maximum allowed chain length

// Uncomment the following macro if you want extra verbose logging in assembly phase.
// #define VERBOSE_ASSEMBLY 1

int main(int argc, char **argv) {
    upcxx::init();

    if (argc < 2) {
        BUtil::print("usage: srun -N nodes -n ranks ./kmer_hash kmer_file [verbose|test [prefix]]\n");
        upcxx::finalize();
        exit(1);
    }

    std::string kmer_fname = std::string(argv[1]);
    std::string run_type = "";
    if (argc >= 3) {
        run_type = std::string(argv[2]);
    }
    std::string test_prefix = "test";
    if (run_type == "test" && argc >= 4) {
        test_prefix = std::string(argv[3]);
    }

    // Check that the kmer file contains kmers of expected length.
    int ks = kmer_size(kmer_fname);
    if (ks != KMER_LEN) {
        throw std::runtime_error("Error: " + kmer_fname + " contains " + std::to_string(ks) +
                                 "-mers, while this binary is compiled for " +
                                 std::to_string(KMER_LEN) + "-mers. Modify packing.hpp and recompile.");
    }

    size_t n_kmers = line_count(kmer_fname);
    size_t hash_table_size = static_cast<size_t>(n_kmers * (1.0 / 0.5));

    upcxx::barrier();
    auto start_insert = std::chrono::high_resolution_clock::now();

    // Create the distributed hash table.
    HashMap hashmap(hash_table_size);

    // Read kmers from the file. Each rank reads its partition.
    std::vector<kmer_pair> kmers = read_kmers(kmer_fname, upcxx::rank_n(), upcxx::rank_me());
    upcxx::barrier();

    // Insertion Phase.
    for (auto& kmer : kmers) {
        bool success = hashmap.insert(kmer);
        if (!success) {
            throw std::runtime_error("Error: HashMap is full during insertion!");
        }
    }
    upcxx::barrier();
    auto end_insert = std::chrono::high_resolution_clock::now();
    double insert_time = std::chrono::duration<double>(end_insert - start_insert).count();

    if (run_type != "test") {
        BUtil::print("Finished inserting in %lf seconds\n", insert_time);
    }

    // Identify starting kmers for contig assembly.
    std::vector<kmer_pair> start_nodes;
    for (const auto &kp : kmers) {
        if (kp.backwardExt() == 'F') { // This condition identifies a contig start.
            start_nodes.push_back(kp);
        }
    }

    // Barrier and time the contig assembly phase.
    upcxx::barrier();
    auto start_assemble = std::chrono::high_resolution_clock::now();

    std::list<std::list<kmer_pair>> contigs;

    // Build contigs from each starting kmer.
    for (const auto& start_kmer : start_nodes) {
        std::list<kmer_pair> contig;
        contig.push_back(start_kmer);

        while (contig.back().forwardExt() != 'F') {
            if (contig.size() >= MAX_CHAIN_LENGTH) {
                BUtil::print("WARNING: Maximum chain length (%d) exceeded for contig starting with %s. Breaking out.\n",
                             MAX_CHAIN_LENGTH, contig.front().kmer.get().c_str());
                break;
            }

            // Check for cycles in the contig.
            pkmer_t next = contig.back().next_kmer();
            bool cycle_found = false;
            for (const auto &existing : contig) {
                if (existing.kmer.get() == next.get()) {
                    cycle_found = true;
                    break;
                }
            }
            if (cycle_found) {
                BUtil::print("WARNING: Cycle detected for contig starting with %s. Breaking out.\n",
                             contig.front().kmer.get().c_str());
                break;
            }

            kmer_pair next_kmer;
            bool found = hashmap.find(next, next_kmer);
            if (!found) {
                throw std::runtime_error("Error: k-mer not found in hashmap during assembly.");
            }
            contig.push_back(next_kmer);
        }
        contigs.push_back(contig);
    }
    upcxx::barrier();
    auto end_assemble = std::chrono::high_resolution_clock::now();
    double assemble_time = std::chrono::duration<double>(end_assemble - start_assemble).count();
    double total_time = std::chrono::duration<double>(end_assemble - start_insert).count();

    int numKmers = std::accumulate(contigs.begin(), contigs.end(), 0,
        [](int sum, const std::list<kmer_pair>& contig) { return sum + contig.size(); });

    if (run_type != "test") {
        BUtil::print("Assembled in %lf seconds total\n", total_time);
        BUtil::print("Constructed %lu contigs with a total of %d nodes from %lu start nodes.\n",
                     contigs.size(), numKmers, start_nodes.size());
    }

    if (run_type == "test") {
        std::ofstream fout(test_prefix + "_" + std::to_string(upcxx::rank_me()) + ".dat");
        for (const auto& contig : contigs) {
            fout << extract_contig(contig) << std::endl;
        }
        fout.close();
    }

    upcxx::finalize();
    return 0;
}