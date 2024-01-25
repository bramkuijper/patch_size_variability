#ifndef _PATCH_SIZE_HPP_
#define _PATCH_SIZE_HPP_

#include <vector>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include "patch.hpp"
#include "parameters.hpp"

class PatchSize
{
    public:
        PatchSize(Parameters const &parms);

    private:
        // random device which is used to generate
        // proper random seeds
        std::random_device rd;
        
        // store the random seed
        // we need to store this so that we can output the
        // random seed, so that we could 'replay' the exact
        // same sequence of random numbers for debugging purposes etc
        unsigned int seed;

        // random number generator
        std::mt19937 rng_r;

        // uniform distribution to compare against probabilities
        std::uniform_real_distribution<double> uniform;

        // uniform distribution to get random patch
        std::uniform_int_distribution<int> patch_sampler;

        // a data file containing the results of the simulation
        std::ofstream data_file;

        // keep track of the time step of the simulation
        long unsigned time_step = 0;

        // parameter object
        // containing all the parameters for this run
        Parameters par;

        // frequency of large patches
        double p_large{0.5};

        // metapopulation of patches
        std::vector<Patch> metapop;

        void initialize_patches();
        void write_parameters();
        void write_data();
        void write_data_headers();

        double mortality(double const zfoc
                , double const zpatch) const;

        void express_help();

        void fill_vacancies();

        // survive or replace vacancy with newborn
        void survive();

        // change patch size
        void change_size();

        // replace breeders 
        void replace();
}; // end class PatchSize

#endif
