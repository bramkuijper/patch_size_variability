#include "patch_size.hpp"
#include "parameters.hpp"

int main(int argc, char **argv)
{
    Parameters params;
    params.d = std::stod(argv[1]);
    params.n[small] = std::stoi(argv[2]);
    params.n[large] = std::stoi(argv[3]);
    params.baseline_mortality = std::stod(argv[4]);
    params.Bfec = std::stod(argv[5]);
    params.Cfec = std::stod(argv[6]);
    params.base_name = argv[7];

    PatchSize sim(params);
}
