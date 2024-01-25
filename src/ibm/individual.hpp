#ifndef INDIVIDUAL_H_
#define INDIVIDUAL_H_

#include <random>
#include "parameters.hpp"

class Individual
{
    public:
        double z[2]{0.0,0.0}; // help in small vs big patches

        Individual(Parameters const &par);

        Individual(Individual const &other);

        // birth constructor
        Individual(Individual const &parent,
                std::mt19937 &rng_r,
                Parameters const &par);
};

#endif
