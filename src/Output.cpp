#include "Output.hpp"

void Tumour::print_to_screen(const InputParameters& params, const DerivedParameters& d_params, const EventCounter& event_counter) {

    std::cout << "Seed = " << params.seed << "; K = " << d_params.K << std::endl;
    std::cout << "Initial migration rate = " << params.init_migration_rate << std::endl;
    std::cout << "Mut rates: " << params.mu_driver_birth << " (birth), " << params.mu_driver_migration << " (mig), " << std::endl;
    std::cout << "Methylation rates: :" << params.meth_rate << " (meth), " << params.demeth_rate << " (demeth)" << std::endl;
    std::cout << "Mut effects: " << params.s_driver_birth << " (birth), " << params.s_driver_migration << " (mig)" << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;

    std::cout << gens_elapsed << " generations, " << iterations << " iterations" << std::endl;
    std::cout << num_cells() << " cells, " << num_clones() << " clones, " << num_driver_genotypes() << " driver genotypes, " << std::endl;
    //std::cout << num_matrix_cols << " matrix columns; " << next_genotype_id << " genotypes ever created" << std::endl;
    std::cout << num_demes() << " demes " <</* max_layer_needed_array[num_demes] + 1 << " bintree layers" <<*/ std::endl;
    std::cout << std::endl;
    std::cout << "----------------------------------------------------------------" << std::endl;

    /*
    std::cout << "Mean number of methylations, demethylations, drivers = " << mean_num_meth << ", " << mean_num_demeth << ", " << mean_num_drivers << std::endl;
    if(calculate_total_diversity) std::cout << "Diversity = " << diversity << " (alpha = " << alpha_diversity << ", beta = " << diversity / alpha_diversity << ")" << std::endl;
    std::cout << "Driver diversity = " << driver_diversity << " (alpha = " << alpha_driver_diversity << ", beta = " << driver_diversity / alpha_driver_diversity << ")" << std::endl;
    if(migration_type < 2) std::cout << "Mean birth, death, mig rates = " << sum_birth_rates / num_cells << ", " << sum_death_rates / num_cells << ", " << sum_migration_rates / num_cells << std::endl;
    else std::cout << "Mean birth, death, fission rates = " << sum_birth_rates / num_cells << ", " << sum_death_rates / num_cells << ", " << sum_migration_rates / num_cells << std::endl;
    if(normal_pop > 0) std::cout << "Mean normal birth, death rates = " << sum_normal_birth_rates / normal_pop << ", " << sum_normal_death_rates / normal_pop << std::endl;
    else std::cout << "No normal cells in tumour demes" << std::endl;

    if(iterations > 0) {
        std::cout << "----------------------------------------------------------------" << std::endl;
        std::cout << event_counter[BIRTH_EVENT] << " births, " << event_counter[DEATH_EVENT] << " deaths, ";
        std::cout << event_counter[FISSION_EVENT] << " fissions, ";
        std::cout << event_counter[MUTATION_EVENT] << " mutations" << std::endl;
        std::cout << event_counter[NORMALBIRTH_EVENT] << " normal cell births, " << event_counter[NORMALDEATH_EVENT] << " normal cell deaths" << std::endl;
    }
    */
}

/*void Tumour::final_output() {
    std::cout << "Placeholder." << std::endl;
}*/