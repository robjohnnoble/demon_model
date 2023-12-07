#include "Objects.hpp"

// Constructor
Clone::Clone(int population, int deme, int genotype, int driverGenotype, int indexInDeme)
    : population(population), deme(deme), genotype(genotype), driver_genotype(driverGenotype), index_in_deme(indexInDeme) {}

// Methods
void Clone::initial_array(const InputParameters& params, const DerivedParameters& d_params, RandomNumberGenerator& rng) {
    if(params.manual_array == -1) {
        for (int i = 0; i < d_params.fcpgs; i++) {
            rng.unif_ran() > 0.5 ? meth_array.push_back(1) : meth_array.push_back(0);
        }
    }
    else {
        for (int i = 0; i < d_params.fcpgs * params.manual_array; i++) {
            meth_array.push_back(0);
        }
        for (int i = std::ceil(d_params.fcpgs * params.manual_array); i < d_params.fcpgs; i++) {
            rng.unif_ran() > 0.5 ? meth_array.push_back(1) : meth_array.push_back(0);
        }
    }
}