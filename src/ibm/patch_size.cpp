#include <stdexcept>
#include <random>
#include <cassert>
#include <vector>
#include <array>
#include "parameters.hpp"
#include "patch_size.hpp"

// constructor
PatchSize::PatchSize(Parameters const &params) :
    rd{} // initialize random device, see *.h file
    ,seed{rd()} // initialize seed
    ,rng_r{seed} // initialize the random number generator
    ,uniform{0.0,1.0} // initialize the uniform distribution
    ,patch_sampler{0,(int)params.npatches - 1} // initialize uniform distribution to sample patch indices from, as we are counting from 0 this can be any number from 0 to npatches - 1
    ,data_file{params.base_name.c_str()} // initialize the data file by giving it a name
    ,par{params} // initialize the parameter data member with the constructor argument
{
    initialize_patches();
    write_data_headers();
    
    for (time_step = 0; time_step <= par.max_time_steps; ++time_step)
    {
        reproduce();

        survive_otherwise_replace();

        if (time_step % par.data_interval == 0 ||
                time_step == par.max_time_steps)
        {
            write_data();
        }
    }
} // end PatchSize()

PatchSize::initialize_patches()
{
    double p_large{par.switch_rate[1] / (par.switch_rate[0] + par.switch_rate[1])};

    for (unsigned patch_idx{0}; patch_idx < par.npatches; ++patch_idx)
    {
        metapop.push_back(Patch(uniform(rng_r) < p_large ? large : small, par));
    }
} // end initialize_patches()


