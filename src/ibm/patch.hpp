#ifndef _PATCH_H_
#define _PATCH_H_

#include <random>
#include "individual.hpp"

enum PatchState
{
    small = 0,
    large = 1
};


class Patch
{
    public:

        PatchState patch_type{small};
        Parameters par{}; 
        double patch_level_z{0.0};

        // two dimensional list of individuals in each patch
        // which reflect the breeders of species 1 and species 2
        std::vector <Individual> breeders;
        std::vector <Individual> breeders_tplus1;
        std::vector <Individual> juveniles;

        // the default constructor 
        // with n1 species 1 breeders
        // with n2 species 2 breeders
        Patch(PatchState const type, Parameters const &params);

        // copy constructor
        Patch(Patch const &other);

        
        void operator=(Patch const &other);

};


#endif
