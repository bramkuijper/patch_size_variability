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
    ,p_large{par.switch_rate[small] / (par.switch_rate[small] + par.switch_rate[large])}
{
    initialize_patches();
    write_data_headers();
    
    for (time_step = 0; time_step <= par.max_time_steps; ++time_step)
    {
        express_help();

        survive();

        change_size();

        fill_vacancies();

        replace();

        if (time_step % par.data_interval == 0 ||
                time_step == par.max_time_steps)
        {
            write_data();
        }
    }

    write_parameters();
} // end PatchSize()

// initialize the metapopuolation 
void PatchSize::initialize_patches()
{
    // go through all patches and initialize
    for (unsigned patch_idx{0}; patch_idx < par.npatches; ++patch_idx)
    {
        metapop.push_back(Patch(uniform(rng_r) < p_large ? large : small, par));
    }
} // end initialize_patches()

void PatchSize::express_help()
{
    double z_patch, z_ind;
    PatchState patch_type;

    // go over all patches and calculate values of z
    for (auto patch_iterator = metapop.begin(); 
            patch_iterator != metapop.end(); 
            ++patch_iterator)
    {
        z_patch = 0;
        patch_type = patch_iterator->patch_type;

        // go through all local breeders and express z
        for (auto ind_iterator = patch_iterator->breeders.begin();
                ind_iterator != patch_iterator->breeders.end();
                ++ind_iterator)
        {
            z_ind = ind_iterator->z[patch_type];

            assert(z_ind >= par.zbound[0]);
            assert(z_ind <= par.zbound[1]);
            z_patch += z_ind;
        }

        z_patch /= patch_iterator->breeders.size();

        patch_iterator->patch_level_z = z_patch;
    }
} // end express_help()

void PatchSize::fill_vacancies()
{
    PatchState patch_type, remote_patch_type;

    std::vector <double> fecundity_weightings_local;
    std::vector <double> fecundity_weightings_remote;

    std::vector < std::pair<unsigned, unsigned> > remote_breeders;

    // whole bunch of aux variables
    double sum_local, sum_remote, fecundity, local_weighting, p_local;

    unsigned n_new_breeders, remote_patch_id, remote_breeder_id;

    // let each patch produce a fixed number of K offspring
    for (auto patch_iterator = metapop.begin(); 
            patch_iterator != metapop.end(); 
            ++patch_iterator)
    {
        // for each newborn build a sampling distribution
        patch_type = patch_iterator->patch_type;

        assert(patch_iterator->breeders_tplus1.size() <= 
                par.n[patch_type]);

        // calculate new offspring
        n_new_breeders = par.n[patch_type] - 
            patch_iterator->breeders_tplus1.size();

        // bounds checking
        assert(n_new_breeders <= par.n[patch_type]);

        fecundity_weightings_local.clear();

        // new recruits needed!
        if (n_new_breeders > 0)
        {
            sum_local = 0.0;

            local_weighting = (double) patch_iterator->breeders.size() / par.n_sample_remote;

            // first get local fecundities
            for (auto individual_iterator = patch_iterator->breeders.begin();
                    individual_iterator != patch_iterator->breeders.end();
                    ++individual_iterator)
            {
                fecundity = local_weighting *
                    (1.0 - par.d) * (1.0 - par.Cfec * individual_iterator->z[patch_type]
                            + par.Bfec * patch_iterator->patch_level_z);

                sum_local += fecundity;

                fecundity_weightings_local.push_back(fecundity);
            }

            std::discrete_distribution<unsigned> local_breeder_sampler(
                    fecundity_weightings_local.begin(),
                    fecundity_weightings_local.end());

            // now start calculating remote fecundities for each newborn
            for (unsigned newborn_idx{0}; 
                    newborn_idx < n_new_breeders; ++newborn_idx)
            {
                // clear the current fecundity weightings
                fecundity_weightings_remote.clear();
                remote_breeders.clear();

                sum_remote = 0.0;

                // sample other breeders from different patches
                for (unsigned remote_idx {0}; 
                        remote_idx < par.n_sample_remote; ++remote_idx)
                {
                    remote_patch_id = patch_sampler(rng_r);

                    std::uniform_int_distribution<unsigned> remote_breeder_sampler(
                            0, metapop[remote_patch_id].breeders.size() - 1);

                    remote_breeder_id = remote_breeder_sampler(rng_r);

                    remote_patch_type = metapop[remote_patch_id].patch_type;

                    // we need to weigh this individual with the frequency
                    // that it would be encountered, which is p_patch_type * n 
                    fecundity = par.d * (double) par.n[remote_patch_type] /
                        ((1.0 - p_large) * par.n[small] + p_large * par.n[large])
                            * (1.0 - metapop[remote_patch_id].breeders[
                                remote_breeder_id
                                ].z[remote_patch_type] * par.Cfec // C * zfoc
                                + 
                                metapop[remote_patch_id].patch_level_z * par.Bfec); 

                    sum_remote += fecundity;

                    fecundity_weightings_remote.push_back(fecundity);

                    remote_breeders.push_back(
                            std::pair<unsigned, unsigned>(remote_patch_id, remote_breeder_id));
                } // end for unsigned remote_idx

                // ok we are done, now make decisions
                p_local = sum_local / (sum_local + sum_remote);

                if (uniform(rng_r) < p_local)
                {
                    Individual Kid(
                            patch_iterator->breeders[local_breeder_sampler(rng_r)],
                            rng_r,
                            par);

                    assert(Kid.z[0] >= par.zbound[0]);
                    assert(Kid.z[0] <= par.zbound[1]);
                    assert(Kid.z[1] >= par.zbound[0]);
                    assert(Kid.z[1] <= par.zbound[1]);

                    patch_iterator->breeders_tplus1.push_back(Kid);
                }
                else
                {
                    // ok more work: we need to sample from remote breeders
                    std::discrete_distribution<unsigned> remote_breeder_sampler{
                        fecundity_weightings_remote.begin(),
                        fecundity_weightings_remote.end()};

                    unsigned remote_idx{remote_breeder_sampler(rng_r)};

                    unsigned remote_patch_idx{
                        std::get<0>(remote_breeders[remote_idx])};

                    unsigned remote_breeder_idx{
                        std::get<1>(remote_breeders[remote_idx])};

                    Individual Kid(
                            metapop[remote_patch_idx].breeders[remote_breeder_idx],
                            rng_r,
                            par);
                    
                    assert(Kid.z[0] >= par.zbound[0]);
                    assert(Kid.z[0] <= par.zbound[1]);
                    assert(Kid.z[1] >= par.zbound[0]);
                    assert(Kid.z[1] <= par.zbound[1]);

                    patch_iterator->breeders_tplus1.push_back(Kid);
                }
            } // end for newborn idx
        } // end if n_new_breeders

        // check whether we have the correct size here
        assert(patch_iterator->breeders_tplus1.size() == 
                par.n[patch_iterator->patch_type]);

    } // end for patch_iterator
} // end fill_vacancies()

void PatchSize::replace()
{
    for (auto patch_iterator = metapop.begin(); 
            patch_iterator != metapop.end(); 
            ++patch_iterator)
    {
        assert(patch_iterator->breeders_tplus1.size() > 0);
        assert(patch_iterator->breeders_tplus1.size() == par.n[patch_iterator->patch_type]);

        patch_iterator->breeders = patch_iterator->breeders_tplus1;
    }
} // end replace()

//  change the size of the patch
void PatchSize::change_size()
{
    // aux variable storing the current type (small, large)
    PatchState current_patch_type, new_patch_type;

    // aux variable storing current count of breeders
    unsigned n_surviving_breeders;

    for (unsigned patch_idx{0}; 
            patch_idx < par.npatches; ++patch_idx)
    {
        current_patch_type = metapop[patch_idx].patch_type;
            
        n_surviving_breeders = 
            metapop[patch_idx].breeders_tplus1.size();

        assert(n_surviving_breeders >= 0);
        assert(n_surviving_breeders <= par.n[current_patch_type]);
            

        // change in size
        if (uniform(rng_r) < par.switch_rate[current_patch_type])
        {
            // switch patch state
            current_patch_type = current_patch_type == small ? large : small;

            // check whether we need to further reduce number of adult 
            // breeders
            if (n_surviving_breeders > par.n[current_patch_type])
            {
                // yes, we do, randomly sample elements to delete
                unsigned n_delete = n_surviving_breeders - par.n[current_patch_type];

                assert(n_delete <= n_surviving_breeders);

                for (unsigned delete_idx{0}; delete_idx < n_delete; ++delete_idx)
                {
                    // random sample breeder to erase
                    std::uniform_int_distribution<size_t> breeder_to_delete{0
                        ,metapop[patch_idx].breeders_tplus1.size() - 1};

                    metapop[patch_idx].breeders_tplus1.erase(
                            metapop[patch_idx].breeders_tplus1.begin()
                            + breeder_to_delete(rng_r));
                }

                assert(metapop[patch_idx].breeders_tplus1.size() 
                        == par.n[current_patch_type]);
            }

            metapop[patch_idx].patch_type = current_patch_type;
        }
    } // end for
} // end change_size()

void PatchSize::survive()
{
    double pmort, patch_level_z;

    PatchState current_patch_state;

    for (auto patch_iter = metapop.begin(); 
            patch_iter < metapop.end(); ++patch_iter)
    {
        current_patch_state = patch_iter->patch_type;
        patch_level_z = patch_iter->patch_level_z;

        // clear existing breeders in breeders_tplus1
        patch_iter->breeders_tplus1.clear();

        // go through all breeders and have them survive
        for (auto it = patch_iter->breeders.begin(); 
                it != patch_iter->breeders.end(); ++it)
        {
            pmort = mortality(
                    it->z[current_patch_state], 
                    patch_level_z);

            // individual will survive
            if (uniform(rng_r) < 1.0 - pmort)
            {
                patch_iter->breeders_tplus1.push_back(*it);
            }
        }
    } // end for patch_idx
} // end survive()

void PatchSize::write_parameters() 
{
    data_file << std::endl << std::endl
        << "z_init;" << par.z_init << std::endl
        << "mu_z0;" << par.mu_z[0] << std::endl
        << "mu_z1;" << par.mu_z[1] << std::endl
        << "seed;" << seed << std::endl
        << "d;" << par.d << std::endl
        << "sdmu;" << par.sdmu << std::endl
        << "n_small;" << par.n[small] << std::endl
        << "n_large;" << par.n[large] << std::endl
        << "Csurv;" << par.Csurv << std::endl
        << "Bsurv;" << par.Bsurv << std::endl
        << "Cfec;" << par.Cfec << std::endl
        << "Bfec;" << par.Bfec << std::endl
        << "baseline_mortality;" << par.baseline_mortality << std::endl
        << "s_small;" << par.switch_rate[small] << std::endl
        << "s_large;" << par.switch_rate[large] << std::endl
        << "p_large;" << p_large << std::endl
        << "n_patches;" << par.npatches << std::endl;
}

void PatchSize::write_data()  
{
    double meanz[2]{0.0,0.0};
    double varz[2]{0.0,0.0};

    double z;

    unsigned n{0};

    unsigned n_large{0};

    // collect stats
    for (auto patch_iterator = metapop.begin(); 
            patch_iterator != metapop.end(); 
            ++patch_iterator)
    {
        for (auto it = patch_iterator->breeders.begin(); 
                it != patch_iterator->breeders.end(); ++it)
        {
            for (unsigned type_idx{0}; type_idx < 2; ++type_idx)
            {
                z = it->z[type_idx];
                meanz[type_idx] += z;
                varz[type_idx] += z * z;
            }
        }

        n += patch_iterator->breeders.size();

        if (patch_iterator->patch_type == large)
        {
            n_large += patch_iterator->breeders.size();
        }
    }

    data_file << time_step << ";";

    // write out the stats
    for (unsigned type_idx{0}; type_idx < 2; ++type_idx)
    {
        meanz[type_idx]/=n;
        varz[type_idx] = varz[type_idx]/n - meanz[type_idx] * meanz[type_idx];

        data_file << (meanz[type_idx]/n) << ";" << varz[type_idx] << ";";
    }

    data_file << n << ";" << (double)n_large / n << std::endl;
}

void PatchSize::write_data_headers() 
{
    data_file << "time;";
    std::string labels[2]{"",""};
    labels[small]="small";
    labels[large]="large";
    
    for (unsigned type_idx{0}; type_idx < 2; ++type_idx)
    {
        data_file << "z_" << labels[type_idx] << ";" << "varz_" << labels[type_idx] << ";";
    }

    data_file << "n;freq_large" << std::endl;
}

// mortality probability
double PatchSize::mortality(double const zfocal, double const zpatch) const
{
    return(par.baseline_mortality
            + (1.0 - par.baseline_mortality) * (
                par.Csurv * zfocal - par.Bsurv * zpatch));
}
