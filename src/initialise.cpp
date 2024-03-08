#include "initialise.hpp"

// set number of generations for the run so that the size of the tumour is
// roughly constant (around 5,000,000 glands)
int calculateMaxGenerations(float fission_rate, int deme_carrying_capacity) {
    return 0;
}

DerivedParameters deriveParameters(const InputParameters& params) {
    DerivedParameters d_params;
    d_params.fcpgs = params.fCpG_loci_per_cell * 2;
    d_params.dmax = 10;
    d_params.max_clones_per_deme = std::ceil(max(d_params.K + 5, 1.2 * d_params.K));
    int max_driver_mu = max(params.mu_driver_migration, params.mu_driver_birth);
    int max_pop = max(9 * d_params.K, 8 * d_params.K + 20);
    d_params.max_genotypes = min(400 * max(max_driver_mu * max_pop, (params.meth_rate + params.demeth_rate) * max_pop), 1e8);
	d_params.max_genotypes = max(d_params.max_genotypes, 10);
    d_params.max_driver_genotypes = min(400 * max(max_driver_mu * max_pop, (params.meth_rate + params.demeth_rate) * max_pop), 1e8);
    int predicted_clones_per_deme = max(min(std::ceil(d_params.K * max_driver_mu * 400), d_params.K), 1);
    d_params.max_clones = min(max(d_params.max_genotypes, 8 * predicted_clones_per_deme), max_pop);
    d_params.max_generations = calculateMaxGenerations(params.init_migration_rate, params.deme_carrying_capacity);
    return d_params;
}
