// compare_files.cpp
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>

int main(int argc, char* argv[]) {
    if(argc < 3) {
        std::cerr << "Usage: " << argv[0] << " file1 file2" << std::endl;
        return 1;
    }

    std::ifstream fin1(argv[1], std::ios::binary);
    std::ifstream fin2(argv[2], std::ios::binary);

    if (!fin1) {
        std::cerr << "Error opening file: " << argv[1] << std::endl;
        return 2;
    }
    if (!fin2) {
        std::cerr << "Error opening file: " << argv[2] << std::endl;
        return 3;
    }

    // Read entire files into strings.
    std::string file1((std::istreambuf_iterator<char>(fin1)), std::istreambuf_iterator<char>());
    std::string file2((std::istreambuf_iterator<char>(fin2)), std::istreambuf_iterator<char>());
    fin1.close();
    fin2.close();

    size_t len1 = file1.size();
    size_t len2 = file2.size();
    size_t min_len = std::min(len1, len2);

    size_t mismatch_count = 0;
    size_t first_mismatch = min_len; // if remains min_len, then no mismatch in common portion

    // Compare character by character.
    for(size_t i = 0; i < min_len; i++){
        if(file1[i] != file2[i]){
            if(first_mismatch == min_len) {
                first_mismatch = i;
            }
            mismatch_count++;
        }
    }

    // Account for extra characters if file lengths differ.
    if(len1 != len2) {
        std::cout << "Files differ in length:" << std::endl;
        std::cout << "  " << argv[1] << " length = " << len1 << std::endl;
        std::cout << "  " << argv[2] << " length = " << len2 << std::endl;
        mismatch_count += std::abs(static_cast<long long>(len1) - static_cast<long long>(len2));
        // If there was no mismatch in the common portion,
        // we set first_mismatch to min_len.
        if(first_mismatch == min_len)
            first_mismatch = min_len;
    }
    
    if(mismatch_count == 0) {
        std::cout << "Files match exactly." << std::endl;
    }
    else {
        std::cout << "Total mismatches: " << mismatch_count 
                  << "\nFirst mismatch at position: " << first_mismatch << std::endl;
    }

    return 0;
}