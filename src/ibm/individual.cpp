#include <random>
#include <algorithm>
#include <set>
#include "individual.hpp"
#include "parameters.hpp"

// default constructor
Individual::Individual(Parameters const &par) :
    z{par.z_init,z_init}
{
} // end Individual()

Individual::Individual(Individual const &parent,
                std::mt19937 &rng_r,
                Parameters const &par) :
    z{parent.z[0], parent.z[1]}
{
    std::uniform_real_distribution uniform{0.0,1.0};

    if (uniform(rng_r) < par.mu_z)
    {
        std::normal_distribution mutational_effect_size{0.0,par.sdmu};
        z = std::clamp(z + mutational_effect_size(rng_r),0.0,1.0);
    }
} // end birth constructor
  
