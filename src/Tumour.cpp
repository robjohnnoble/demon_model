#include "Tumour.hpp"

// initialise tumour from input parameters
void Tumour::initialise(const InputParameters& params,
    const DerivedParameters& d_params, RandomNumberGenerator& rng) {
    // genotypes:
    Genotype genotype(params.init_pop, 0, 0, 0, 0, 0, 0, 0, 1, 1, params.init_migration_rate, 0);
    genotypes.push_back(genotype);

    // driver genotypes:
    DriverGenotype driver_genotype(params.init_pop, 0, 0, 0, 0, 0, 0, 0, 1, 1, params.init_migration_rate, 0);
    driver_genotypes.push_back(driver_genotype);

    // clones:
    Clone clone(1, 0, 0, 0, 0);
    clone.initial_array(params, d_params, rng);
    clones.push_back(clone);

    // demes:
    Deme deme(d_params.K, "left", 0, 1, 1, 0, 0, 0, 1, params.init_migration_rate);
    deme.initial_sum();
    deme.calculate_average_array(d_params, clones);

    // fission times
    fission_times.push_back(params.time0);
    fission_times.push_back(params.time1);
    fission_times.push_back(params.time2);
    fission_times.push_back(params.time3);
    fission_times.push_back(params.time4);
    fission_times.push_back(params.time5);
    fission_times.push_back(params.time6);

    next_fission = fission_times[0];
}

// deme fission
void Tumour::deme_fission(){
    if (gens_elapsed >= next_fission) {
        // fission occurs
    }
    else {
        // pseudo fission
    }
}

// division of cancer cells
void Tumour::cell_division(int parent_clone, EventCounter& event_counter, RandomNumberGenerator& rng,
    Deme& deme, Clone& clone, Genotype& genotype, DriverGenotype& driver_genotype, const InputParameters& params) {
    int daughter_clone_nums[2] = {parent_clone, parent_clone};

    event_counter.birth++;
    deme.increment(1, params);

    std::vector<int> new_birth_drivers(2,0), new_mig_drivers(2,0);
    std::vector<int> new_mutations = choose_number_mutations(rng, params.mu_driver_birth, params.mu_driver_migration, new_birth_drivers, new_mig_drivers);

    if(!new_mutations[0] && !new_mutations[1]) {
        driver_genotype.increment(1);
    }

}

// create new genotype upon division
void Tumour::create_genotype(Genotype& parent) {
    Genotype genotype(1, parent.identity, next_genotype_id++, parent.driver_identity, parent.number_of_driver_mutations, parent.number_of_migration_mutations,
    parent.num_meth, parent.num_demeth, 1, parent.birth_rate, parent.migration_rate, gens_elapsed);
}

// mutation event handling
std::vector<int> Tumour::choose_number_mutations(RandomNumberGenerator& rng, float mu_driver_birth, float mu_driver_migration, 
    std::vector<int>& new_birth_drivers, std::vector<int>& new_mig_drivers) {
    std::vector<int> res(2, 0);
    for(int i = 0; i < 2; i++) {
        new_birth_drivers[i] = rng.poisson(mu_driver_birth);
        new_mig_drivers[i] = rng.poisson(mu_driver_migration);
        res[i] = new_birth_drivers[i] + new_mig_drivers[i];
    }
    return res;
}

// calculate sums of different rates in the whole tumour
void Tumour::calculate_sums_of_rates() {
    for (int i = 0; i < demes.size(); i++) {
        sum_death_rates += demes[i].death_rate * demes[i].population;
        sum_birth_rates += demes[i].sum_birth_rates;
        sum_migration_rates += demes[i].sum_migration_rates;
    }
}

double Tumour::sum_of_all_rates() {
    return sum_birth_rates + sum_death_rates + sum_migration_rates;
}

// update time
void Tumour::update_time(RandomNumberGenerator& rng, const InputParameters& params) {
    float temp_sum = sum_of_all_rates();

    float gens_added = rng.exp(1);
    gens_elapsed += gens_added;
    output_timer += gens_added;
}

//check time
float Tumour::check_time() {
    return gens_elapsed;
}

// numbers of relevant things
int Tumour::num_cells() {
    int res = 0;
    for (int i = 0; i < demes.size(); i++) {
        res += demes[i].population;
    }
    return res;
}

int Tumour::num_clones() {
    return clones.size();
}

int Tumour::num_driver_genotypes() {
    return driver_genotypes.size();
}

int Tumour::num_genotypes() {
    return genotypes.size();
}

int Tumour::num_demes() {
    return demes.size();
}