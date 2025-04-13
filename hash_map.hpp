#pragma once
#include <upcxx/upcxx.hpp>
#include <vector>
#include "kmer_t.hpp"
#include "butil.hpp"

struct HashMap {
    upcxx::global_ptr<kmer_pair> local_data;
    upcxx::global_ptr<int> local_used;

    std::vector<upcxx::global_ptr<kmer_pair>> all_data_ptrs;
    std::vector<upcxx::global_ptr<int>> all_used_ptrs;

    size_t my_size;
    size_t global_size;
    size_t local_size;

    int nranks;
    int me;

    HashMap(size_t size);
    bool insert(const kmer_pair &kmer);
    bool find(const pkmer_t &key, kmer_pair &val);

    size_t size() const noexcept { return global_size; }
};

inline int owner_of_slot(uint64_t slot, int nranks, size_t local_size) {
    return static_cast<int>(slot / local_size);
}
inline uint64_t offset_of_slot(uint64_t slot, size_t local_size) {
    return slot % local_size;
}

HashMap::HashMap(size_t size) {
    upcxx::barrier();

    nranks = upcxx::rank_n();
    me = upcxx::rank_me();

    global_size = size;
    local_size = (global_size + nranks - 1) / nranks;
    my_size = local_size;

    local_data = upcxx::new_array<kmer_pair>(my_size);
    local_used = upcxx::new_array<int>(my_size);

    kmer_pair* dptr = local_data.local();
    int* uptr = local_used.local();
    for (size_t i = 0; i < my_size; i++) {
        dptr[i] = kmer_pair();   // Initialize to default (empty bucket)
        uptr[i] = 0;
    }

    all_data_ptrs.resize(nranks);
    all_used_ptrs.resize(nranks);
    all_data_ptrs[me] = local_data;
    all_used_ptrs[me] = local_used;

    upcxx::barrier();
    for (int r = 0; r < nranks; r++) {
        all_data_ptrs[r] = upcxx::broadcast(all_data_ptrs[r], r).wait();
        all_used_ptrs[r] = upcxx::broadcast(all_used_ptrs[r], r).wait();
    }
    upcxx::barrier();
}

bool HashMap::insert(const kmer_pair &kmer) {
    uint64_t hashv = kmer.hash();
    for (uint64_t probe = 0; probe < global_size; probe++) {
        uint64_t slot = (hashv + probe) % global_size;
        int owner = owner_of_slot(slot, nranks, local_size);
        uint64_t off = offset_of_slot(slot, local_size);

        bool claimed = upcxx::rpc(owner, [off, used_ptr = all_used_ptrs[owner]]() {
            int* used = used_ptr.local();
            if (used[off] == 0) {
                used[off] = 1;
                return true;
            }
            return false;
        }).wait();

        if (claimed) {
            upcxx::rpc(owner, [off, kmer, data_ptr = all_data_ptrs[owner]]() {
                kmer_pair* data = data_ptr.local();
                data[off] = kmer;
            }).wait();
            return true;
        }
    }
    return false;
}

bool HashMap::find(const pkmer_t &key, kmer_pair &val) {
    uint64_t hashv = key.hash();
    for (uint64_t probe = 0; probe < global_size; probe++) {
        uint64_t slot = (hashv + probe) % global_size;
        int owner = owner_of_slot(slot, nranks, local_size);
        uint64_t off = offset_of_slot(slot, local_size);

        bool used = upcxx::rpc(owner, [off, used_ptr = all_used_ptrs[owner]]() {
            return (used_ptr.local()[off] != 0);
        }).wait();

        if (!used)
            return false;

        kmer_pair occupant = upcxx::rpc(owner, [off, data_ptr = all_data_ptrs[owner]]() {
            return data_ptr.local()[off];
        }).wait();

        // Compare using the string obtained via get() for each key.
        if (occupant.kmer.get() == key.get()) {
            val = occupant;
            return true;
        }
    }
    return false;
}