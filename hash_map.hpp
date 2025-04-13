#pragma once
#include <upcxx/upcxx.hpp>
#include <vector>
#include <stdexcept>
#include <string>
#include "kmer_t.hpp"
#include "butil.hpp"

// Utility functions to compute owner and offset from a global slot index.
inline int owner_of_slot(uint64_t slot, int nranks, size_t local_size) {
    return static_cast<int>(slot / local_size);
}
inline uint64_t offset_of_slot(uint64_t slot, size_t local_size) {
    return slot % local_size;
}

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

    // Constructor allocates the local arrays on each rank
    // and broadcasts the pointers so that every rank can access any portion.
    HashMap(size_t size) {
        upcxx::barrier();

        nranks = upcxx::rank_n();
        me = upcxx::rank_me();

        global_size = size;
        local_size = (global_size + nranks - 1) / nranks;
        my_size = local_size;

        local_data = upcxx::new_array<kmer_pair>(my_size);
        local_used = upcxx::new_array<int>(my_size);

        // Initialize local arrays.
        kmer_pair* dptr = local_data.local();
        int* uptr = local_used.local();
        for (size_t i = 0; i < my_size; i++) {
            dptr[i] = kmer_pair();  // Default initialize each bucket
            uptr[i] = 0;
        }

        // Prepare vectors to hold remote pointers.
        all_data_ptrs.resize(nranks);
        all_used_ptrs.resize(nranks);
        all_data_ptrs[me] = local_data;
        all_used_ptrs[me] = local_used;

        upcxx::barrier();
        // Broadcast each rank's pointer so that every rank has the complete directory.
        for (int r = 0; r < nranks; r++) {
            all_data_ptrs[r] = upcxx::broadcast(all_data_ptrs[r], r).wait();
            all_used_ptrs[r] = upcxx::broadcast(all_used_ptrs[r], r).wait();
        }
        upcxx::barrier();
    }

    // Insert the provided kmer_pair into the hash table.
    // Uses linear probing and one‐sided RPC to claim a slot.
    bool insert(const kmer_pair &kmer) {
        uint64_t hashv = kmer.hash();
        // (Optional) Uncomment to print insertion info:
        // BUtil::print("Rank %d: Inserting key %s, hash %llu\n",
        //          upcxx::rank_me(), kmer.kmer.get().c_str(), (unsigned long long)hashv);

        for (uint64_t probe = 0; probe < global_size; probe++) {
            uint64_t slot = (hashv + probe) % global_size;
            int owner = owner_of_slot(slot, nranks, local_size);
            uint64_t off = offset_of_slot(slot, local_size);

            // (Optional) Commented out detailed debug:
            // BUtil::print("Rank %d: Insert probe %llu: slot %llu, owner %d, offset %llu\n",
            //              upcxx::rank_me(), (unsigned long long)probe,
            //              (unsigned long long)slot, owner, (unsigned long long)off);

            bool claimed = upcxx::rpc(owner,
                [off, used_ptr = all_used_ptrs[owner]]() mutable -> bool {
                    int *used = used_ptr.local();
                    if (used[off] == 0) {
                        used[off] = 1;
                        return true;
                    }
                    return false;
                }).wait();

            if (claimed) {
                upcxx::rpc(owner,
                    [off, kmer, data_ptr = all_data_ptrs[owner]]() mutable {
                        kmer_pair* data = data_ptr.local();
                        data[off] = kmer;
                    }).wait();
                // BUtil::print("Rank %d: Insert successful at slot %llu (owner %d, offset %llu)\n",
                //              upcxx::rank_me(), (unsigned long long)slot, owner, (unsigned long long)off);
                return true;
            }
        }
        return false;
    }

    // Find the kmer_pair whose key matches the given pkmer_t.
    // Uses linear probing with one‐sided RPCs to retrieve the stored value.
    bool find(const pkmer_t &key, kmer_pair &val) {
        std::string key_str = key.get();
        uint64_t hashv = key.hash();
        // BUtil::print("Rank %d: Looking up key %s, hash %llu\n",
        //          upcxx::rank_me(), key_str.c_str(), (unsigned long long)hashv);

        for (uint64_t probe = 0; probe < global_size; probe++) {
            uint64_t slot = (hashv + probe) % global_size;
            int owner = owner_of_slot(slot, nranks, local_size);
            uint64_t off = offset_of_slot(slot, local_size);

            // BUtil::print("Rank %d: Find probe %llu: slot %llu, owner %d, offset %llu\n",
            //              upcxx::rank_me(), (unsigned long long)probe,
            //              (unsigned long long)slot, owner, (unsigned long long)off);

            bool used = upcxx::rpc(owner,
                [off, used_ptr = all_used_ptrs[owner]]() -> bool {
                    return (used_ptr.local()[off] != 0);
                }).wait();

            if (!used) {
                // BUtil::print("Rank %d: Slot %llu is empty during lookup for key %s\n",
                //              upcxx::rank_me(), (unsigned long long)slot, key_str.c_str());
                return false;
            }

            kmer_pair occupant = upcxx::rpc(owner,
                [off, data_ptr = all_data_ptrs[owner]]() -> kmer_pair {
                    return data_ptr.local()[off];
                }).wait();

            if (occupant.kmer.get() == key.get()) {
                val = occupant;
                return true;
            }
        }
        return false;
    }
};