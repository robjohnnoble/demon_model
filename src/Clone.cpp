#include "Objects.hpp"

// Constructor
Clone::Clone(int population, int deme, int genotype, int driverGenotype, int indexInDeme, int driver_index, int index)
    : population(population), deme(deme), genotype(genotype), driver_genotype(driverGenotype), index_in_deme(indexInDeme), driver_index(driver_index), index(index) {}

// Methods
void Clone::initial_array(const InputParameters& params, const DerivedParameters& d_params, RandomNumberGenerator& rng) {
    meth_array = std::vector<int>(d_params.fcpgs, 0);
    if(params.manual_array == -1) {
        for (int i = 0; i < d_params.fcpgs; i++) {
            rng.unif_ran() > 0.5 ? meth_array[i] = 1 : meth_array [i] = 0;
        }
    }
    else {
        for (int i = 0; i < d_params.fcpgs * params.manual_array; i++) {
            meth_array[i] = 0;
        }
        for (int i = std::ceil(d_params.fcpgs * params.manual_array); i < d_params.fcpgs; i++) {
            rng.unif_ran() > 0.5 ? meth_array[i] = 1 : meth_array [i] = 0;
        }
    }
}