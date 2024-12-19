#include "../helper_functions.hpp"
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <cstdlib>

int main(int argc, char* argv[]){
    int n = 2;
    if (argc > 1) n = std::atoi(argv[1]);
    int N = n * n;

    // initialising output
    std::string filename = std::to_string(n) + "_energy_density.csv";
    std::ofstream output_file(filename);
    output_file << "E,N(E)\n";

    // initialising important variables
    std::vector<std::array<int, 4>> neighbours = createNeighbourList(n); // neighbour list
    unsigned long long total_config = (1ULL << (N-1)) - 1; // up-down symmetry so N-1
    std::vector<unsigned long long> energy_count(4 * N, 0ULL); // vector indexed by energy + 2N to store count

    unsigned long long gray_p = 0; // initial configuration
    int energy = 0; // initial energy + 2N
    energy_count[0] = 2; // adding ground state

    // main loop
    auto start = std::chrono::high_resolution_clock::now();

    for (unsigned long long p = 0; p < total_config; ++p) {
        // finding the flipped spin, index is equivalent to CTZ((p^(p >> 1)) ^ ((p+1) ^ ((p+1) >> 1))) = CTZ(p+1)
        int flip_bit = __builtin_ctzll(p+1); 
        bool flip_bool = (gray_p >> flip_bit) & 1; // finding the value of the flipped spin
        
        // calculating neigbour interactions
        for (int neighbour : neighbours[flip_bit]) {
            energy += 1 - 2 * (((gray_p >> neighbour) & 1) ^ flip_bool);
        }

        gray_p ^= (1ULL << flip_bit); // flipping the spin
        energy_count[energy] += 2; // because of updown symmetry we add two
    }

    auto end = std::chrono::high_resolution_clock::now();

    // writing to file
    for (int E = 0; E <= 4 * N; E++) {
    unsigned long long count = energy_count[E];
    if (count > 0)
        output_file << 2 * E - 2 * N << "," << count << "\n";
    }
    
    output_file.close();

    // printing time taken
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Processing completed in " << elapsed.count() << " seconds.\n";
    std::cout << "Results written to " << filename << "\n";

    return 0;
}