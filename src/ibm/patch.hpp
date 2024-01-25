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
        // two dimensional list of individuals in each patch
        // which reflect the breeders of species 1 and species 2
        std::vector < Individual > breeders;
        std::vector < Individual > juveniles;

        // the default constructor 
        // with n1 species 1 breeders
        // with n2 species 2 breeders
        Patch(int const n, Parameters const &params);

        // copy constructor
        Patch(Patch const &other);

        PatchState patch_type{small};

        Parameters par; 
        
        void operator=(Patch const &other);

};


#endif
