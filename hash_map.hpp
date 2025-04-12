#pragma once

#include "kmer_t.hpp"
#include <map>
#include <upcxx/upcxx.hpp>

struct HashMap {
    std::vector<kmer_pair> data;
    std::vector<int> used;
    using dist_map = upcxx::dist_object<std::vector<std::vector<kmer_pair>>>;
    dist_map local_map;


    int count = 0;

    size_t global_size;

    size_t size() const noexcept;

    HashMap(size_t size);

    // Most important functions: insert and retrieve
    // k-mers from the hash table.
    bool insert(const kmer_pair& kmer);
    bool find(const pkmer_t& key_kmer, kmer_pair& val_kmer);

    // Helper functions

    bool find_in_bucket(uint64_t local_idx, kmer_pair& val_kmer);

    // Write and read to a logical data slot in the table.
    void write_slot(uint64_t slot, const kmer_pair& kmer);
    kmer_pair read_slot(uint64_t slot);

    // Request a slot or check if it's already used.
    bool request_slot(uint64_t slot);
    bool slot_used(uint64_t slot);
};

HashMap::HashMap(size_t size) : global_size(size), local_map(std::vector<std::vector<kmer_pair>>(size / upcxx::rank_n() + 1)){
    // my_size = size;
    // local_map.resize(size);

    // data.resize(size);
    // used.resize(size, 0);
}

bool HashMap::insert(const kmer_pair& kmer) {
    uint64_t global_idx = kmer.hash() % global_size;
    uint64_t target_rank = kmer.hash() % upcxx::rank_n();
    uint64_t local_idx = global_idx % (global_size / upcxx::rank_n() + 1);

    // std::cout << "Hello world from process " << upcxx::rank_me()
    // << " out of " << upcxx::rank_n() << " processes" << std::endl;

    // std::cout << target_rank << std::endl;

    if (target_rank == upcxx::rank_me()) {
        (*local_map)[local_idx].push_back(kmer);
        ++count;
        // if (upcxx::rank_me() == 0)
        // std::cout << count << " on processor " << upcxx::rank_me() << std::endl;
    } else {
        upcxx::rpc(target_rank, [local_idx, kmer](dist_map &dmap) {
            dmap->at(local_idx).push_back(kmer);
        }, local_map).wait();
        ++count;
        // if (upcxx::rank_me() == 0)
        // std::cout << count << " on processor " << upcxx::rank_me() << std::endl;
    }
    return true;

    // uint64_t hash = kmer.hash();
    // uint64_t probe = 0;
    // bool success = false;
    // do {
    //     uint64_t slot = (hash + probe++) % size();
    //     success = request_slot(slot);
    //     if (success) {
    //         write_slot(slot, kmer);
    //     }
    // } while (!success && probe < size());
    // return success;
}

bool HashMap::find(const pkmer_t& key_kmer, kmer_pair& val_kmer) {
    uint64_t global_idx = key_kmer.hash() % global_size;
    uint64_t target_rank = key_kmer.hash() % upcxx::rank_n();
    uint64_t local_idx = global_idx % (global_size / upcxx::rank_n() + 1);

    if (target_rank == upcxx::rank_me()) {
        std::vector<kmer_pair> bucket = (*local_map)[local_idx];
        // for (int i = 0; i < (*local_map).size(); ++i) {
        //     std::vector<kmer_pair> b = (*local_map)[i];
        //     for (int j = 0; j < b.size(); ++j) {
        //         if (b[j] == val_kmer) {
        //             std::cout << "bucket: " << b[j].kmer_str() << " val: " << val_kmer.kmer_str() << std::endl;
        //         }
        //     }
        // }
        for (int i = 0; i < bucket.size(); ++i) {
            std::cout << "bucket: " << bucket[i].kmer_str() << " val: " << val_kmer.kmer_str() << std::endl;
            std::cout << "bucket hash: " << bucket[i].hash() << " key hash: " << key_kmer.hash() << std::endl;
            if (bucket[i] == val_kmer || bucket.size() == 1) {
                val_kmer = bucket[i];
                std::cout << "FOUND :)" << std::endl;
                return true;
            }
        }
        return false;
        // return find_in_bucket(local_idx, val_kmer);
    } else {
        return upcxx::rpc(target_rank, 
            [local_idx, &val_kmer](dist_map& remote) -> bool {
                std::vector<kmer_pair> bucket = (*remote)[local_idx];
                for (int i = 0; i < bucket.size(); ++i) {
                    std::cout << "iterating remotely" << std::endl;
                    if (bucket[i] == val_kmer || bucket.size() == 1) {
                        val_kmer = bucket[i];
                        std::cout << "FOUND :)" << std::endl;
                        return true;
                    }
                }
                return false;
            }, local_map).wait();
    }

    return true;
    // uint64_t hash = key_kmer.hash();
    // uint64_t probe = 0;
    // bool success = false;
    // do {
    //     uint64_t slot = (hash + probe++) % size();
    //     if (slot_used(slot)) {
    //         val_kmer = read_slot(slot);
    //         if (val_kmer.kmer == key_kmer) {
    //             success = true;
    //         }
    //     }
    // } while (!success && probe < size());
    // return success;
}

bool HashMap::find_in_bucket(uint64_t local_idx, kmer_pair& val_kmer) {
    std::vector<kmer_pair> bucket = (*local_map)[local_idx];
    for (int i = 0; i < bucket.size(); ++i) {
        if (bucket[i] == val_kmer) {
            val_kmer = bucket[i];
            return true;
        }
    }
    return false;
}

bool HashMap::slot_used(uint64_t slot) { return used[slot] != 0; }

void HashMap::write_slot(uint64_t slot, const kmer_pair& kmer) { data[slot] = kmer; }

kmer_pair HashMap::read_slot(uint64_t slot) { return data[slot]; }

// need atomic
bool HashMap::request_slot(uint64_t slot) {
    if (used[slot] != 0) {
        return false;
    } else {
        used[slot] = 1;
        return true;
    }
}

size_t HashMap::size() const noexcept { return global_size; }