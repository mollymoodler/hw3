#pragma once

#include "kmer_t.hpp"
#include <upcxx/upcxx.hpp>

struct HashMap {
    using dist_map = upcxx::dist_object<std::vector<std::vector<kmer_pair>>>;
    dist_map local_map;

    size_t global_size;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);
};

HashMap::HashMap(size_t size) : global_size(size), local_map(std::vector<std::vector<kmer_pair>>(size / upcxx::rank_n() + 1)) {}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t global_idx = kmer.hash() % global_size;
    uint64_t target_rank = kmer.hash() % upcxx::rank_n();
    uint64_t local_idx = global_idx % (global_size / upcxx::rank_n() + 1);

    if (target_rank == upcxx::rank_me()) {
        (*local_map)[local_idx].push_back(kmer);
    } else {
        upcxx::rpc(target_rank, [local_idx, kmer](dist_map &dmap) {
            dmap->at(local_idx).push_back(kmer);
        }, local_map).wait();
    }
    return true;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t global_idx = key_kmer.hash() % global_size;
    uint64_t target_rank = key_kmer.hash() % upcxx::rank_n();
    uint64_t local_idx = global_idx % (global_size / upcxx::rank_n() + 1);

    if (target_rank == upcxx::rank_me()) {
        std::vector<kmer_pair> bucket = (*local_map)[local_idx];
        for (int i = 0; i < bucket.size(); ++i) {
            if (bucket[i].kmer == key_kmer) {
                val_kmer = bucket[i];
                return true;
            }
        }
        return false;
    } else {
        kmer_pair val = upcxx::rpc(target_rank, 
            [target_rank, local_idx, key_kmer](dist_map& remote) -> kmer_pair {
                std::vector<kmer_pair> bucket = (*remote)[local_idx];
                for (int i = 0; i < bucket.size(); ++i) {
                    if (bucket[i].kmer == key_kmer) {
                        return bucket[i];
                    }
                }
                return kmer_pair();
            }, local_map).wait();
        if (val.kmer == key_kmer) {
            val_kmer = val;
            return true;
        }
    }

    return true;
}