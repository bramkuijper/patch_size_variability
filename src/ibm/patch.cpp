#include "patch.h"
#include "individual.h"
#include <random>
#include <vector>
#include <cassert>
#include <algorithm>

Patch::Patch(PatchState const type // whether type small or large
        ,Parameters const &params) : // parameter object
    patch_type{type} // initialize patch type data member
    ,par{params} // initialize parameter object data member
    breeders(par.n[patch_type], Individual(par)) // assign n_i breeders to this patch
{}

// copy constructor
Patch::Patch(Patch const &other) :
    patch_type{other.patch_type} 
    ,par{other.par} 
    ,breeders{other.breeders}
{}

// assignment operator
void Patch::operator=(Patch const &other)
{
    patch_type = other.patch_type;
    par = other.par;
    breeders = other.breeders;
} // end assignment operator
