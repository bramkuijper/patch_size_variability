#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP

#include <string>

class Parameters
{
    public:
        double z_init{0.0};
        double mu_z[2]{0.01,0.01};
        double sdmu{0.01};

        // patch sizes of small and large patches
        unsigned n[2]{2,5};

        // sample remote breeders when calculating fecundity
        unsigned n_sample_remote{10};

        // cost of survival help 
        double Csurv{1.0};

        // benefit coefficient of fecundity help
        double Bsurv{1.0};

        double baseline_mortality{0.05};
        
        // cost of survival help 
        double Cfec{1.0};

        // benefit coefficient of fecundity help
        double Bfec{1.0};

        // dispersal probability
        double d{0.2};

        unsigned npatches{500};
        long unsigned max_time_steps{100000};
        long unsigned data_interval{10};

        double switch_rate[2]{0.2,0.9};

        std::string base_name{"sim_patch_size_vary"};

};

#endif
