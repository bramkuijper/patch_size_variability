#ifndef _PARAMETERS_HPP
#define _PARAMETERS_HPP


class Parameters
{
    public:
        double z_init{0.0};
        double mu_z{0.01};
        double sdmu{0.01};

        // patch sizes of small and large patches
        unsigned n[2]{2,5};

        double switch_rate[2]{0.2,0.9};

};

#endif
