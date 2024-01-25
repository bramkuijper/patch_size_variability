#include <random>
#include <algorithm>
#include <set>
#include "individual.hpp"
#include "parameters.hpp"

// default constructor
Individual::Individual(Parameters const &par) :
    z{par.z_init,par.z_init}
{
} // end Individual()

// copy constructor
Individual::Individual(Individual const &other) :
    z{other.z[0],other.z[1]}
{
} // end Individual()

// birth constructory
Individual::Individual(Individual const &parent,
                std::mt19937 &rng_r,
                Parameters const &par) :
    z{parent.z[0], parent.z[1]}
{
    std::uniform_real_distribution uniform{0.0,1.0};

    for (unsigned locus_idx = 0; locus_idx < 2; ++locus_idx)
    {
        if (uniform(rng_r) < par.mu_z[locus_idx])
        {
            std::normal_distribution mutational_effect_size{0.0,par.sdmu};
            z[locus_idx] = std::clamp(z[locus_idx] + mutational_effect_size(rng_r),0.0,1.0);
        }
    }
} // end birth constructor


// assignment operator
void Individual::operator=(Individual const &other)
{
    z[0] = other.z[0];
    z[1] = other.z[1];
}
