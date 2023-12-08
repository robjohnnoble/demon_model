#include "Tumour.hpp"

// initialise tumour from input parameters
void Tumour::initialise(const InputParameters& params,
    const DerivedParameters& d_params, RandomNumberGenerator& rng) {
    // driver genotypes:
    DriverGenotype driver_genotype(params.init_pop, 0, 0, 0, 0, 0, 0, true, 1, params.init_migration_rate, 0);
    driver_genotypes.push_back(driver_genotype);

    // clones:
    Clone clone(1, 0, 0, 0, 0);
    clone.initial_array(params, d_params, rng);
    clones.push_back(clone);

    // demes:
    Deme deme(d_params.K, "left", 0, 1, 0, params.baseline_death_rate, 0, 1, params.init_migration_rate);
    deme.calculate_sum_of_rates();
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

// choose deme
int Tumour::choose_deme(RandomNumberGenerator& rng) {
    std::vector<double> cum_rates(demes.size());
    double r;

    if (cum_rates.size() == 1) {
        return 0;
    } else {
        cum_rates[0] = demes[0].sum_rates;
        for (int i = 1; i < cum_rates.size(); i++) {
            cum_rates[i] = demes[i].sum_rates + cum_rates[i - 1];
        }
        r = rng.unif_ran() * cum_rates.back();
    }
    if (cum_rates.size() == 2) {
        return r < cum_rates[0] ? 0 : 1;
    } else {
        return std::lower_bound(cum_rates.begin(), cum_rates.end(), r) - cum_rates.begin();
    }
}
// choose clone
int Tumour::choose_clone(int chosen_deme, RandomNumberGenerator& rng) {
    // cumulative sum of birth and death rates to choose a clone
    // first need to choose a genotype then a clone which has that driver genotype
    std::vector<double> cum_rates;
    double r;

    if (demes[chosen_deme].clones_list.size() == 1) {
        return 0;
    } else {
        int driver_id = demes[chosen_deme].clones_list[0]->driver_index;
        cum_rates[0] = demes[chosen_deme].death_rate + driver_genotypes[driver_id].birth_rate;
        for (int i = 1; i < demes[chosen_deme].clones_list.size(); i++) {
            cum_rates[i] = cum_rates + demes[chosen_deme].death_rate + driver_genotypes[driver_id].birth_rate;
        }
        r = rng.unif_ran() * cum_rates.back();
    }

    if (demes[chosen_deme].clones_list.size() == 2) {
        return r < cum_rates[0] ? 0 : 1;
    } else {
        return std::lower_bound(cum_rates.begin(), cum_rates.end(), r) - cum_rates.begin();
    }
}

// deme fission
void Tumour::deme_fission(Deme& deme, EventCounter& event_counter, RandomNumberGenerator& rng, const InputParameters& params){
    if (next_fission == fission_times[0]) {
        // first fission - sets up left and right sides of the tumour
        event_counter.fission++;
        Deme new_deme(deme.K, std::string("right"), static_cast<int>(demes.size()), 0, deme.fissions + 1, 0, 0, 0, 0);
        move_cells(deme, new_deme, rng, params);
    } else if (gens_elapsed >= next_fission) {
        // subsequent fissions
        event_counter.fission++;
        Deme new_deme(deme.K, deme.side, static_cast<int>(demes.size()), 0, deme.fissions + 1, 0, 0, 0, 0);
        move_cells(deme, new_deme, rng, params);
    } else {
        // pseudo fission - kills half the population in a deme without creating a new deme
        event_counter.fission++;
        deme.fissions++;
        pseudo_fission(deme, rng, params);
    }
}
// move cells after deme fission
void Tumour::move_cells(Deme& parent, Deme& daughter, RandomNumberGenerator& rng, const InputParameters& params) {
    int cells_to_move = rng.stochastic_round(parent.population / 2);
    
    // create and shuffle a list of indices
    std::vector<int> indices(parent.population);
    std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, 2, ..., parent.population - 1
    std:: shuffle(indices.begin(), indices.end(), rng.get_engine());
    
    // select the first cells_to_move indices and sort by descending order to avoid index errors
    std::sort(indices.begin(), indices.begin() + cells_to_move, std::greater<int>());

    // move the corresponding clones
    for (int i = 0; i < cells_to_move; i++) {
        int index = indices[i];
        daughter.clones_list.push_back(parent.clones_list[index]);
        daughter.clones_list.back()->deme = daughter.identity;
        daughter.clones_list.back()->index_in_deme = daughter.clones_list.size() - 1;
        clones[daughter.clones_list.back()->index].deme = daughter.identity;
        parent.clones_list.erase(parent.clones_list.begin() + index);
    }
    parent.increment(-cells_to_move, params);
    daughter.increment(cells_to_move, params);
}
// pseudo fission
void Tumour::pseudo_fission(Deme& parent, RandomNumberGenerator& rng, const InputParameters& params) {
    int cells_to_move = rng.stochastic_round(parent.population / 2);
    
    // create and shuffle a list of indices
    std::vector<int> indices(parent.population);
    std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, 2, ..., parent.population - 1
    std::shuffle(indices.begin(), indices.end(), rng.get_engine());
    
    // select the first cells_to_move indices and sort by descending order to avoid index errors
    std::sort(indices.begin(), indices.begin() + cells_to_move, std::greater<int>());

    // move the corresponding clones
    for (int i = 0; i < cells_to_move; i++) {
        int index = indices[i];
        remove_clone(parent, *parent.clones_list[index]);
    }
    parent.increment(-cells_to_move, params);
    calculate_deme_birth_rate(parent);
    calculate_deme_migration_rate(parent);
    parent.set_death_rate(params);
}

// remove clone from deme and clones list
void Tumour::remove_clone(Deme& deme, Clone& clone) {
    // remove clone from deme
    deme.clones_list.erase(deme.clones_list.begin() + clone.index_in_deme);
    // remove clone from clones
    clones.erase(clones.begin() + clone.index);
    // update indices
    for (int i = 0; i < deme.clones_list.size(); i++) {
        clones[i].index_in_deme = i;
    }
    for (int i = 0; i < clones.size(); i++) {
        clones[i].index = i;
    }
}
// remove driver genotype from driver genotypes
void Tumour::remove_driver_genotype(DriverGenotype& driver_genotype) {
    driver_genotypes.erase(driver_genotypes.begin() + driver_genotype.index);
    // update indices
    for (int i = 0; i < driver_genotypes.size(); i++) {
        driver_genotypes[i].index = i;
        clones[i].driver_index = i;
    }
}


// division of cancer cells
void Tumour::cell_division(EventCounter& event_counter, RandomNumberGenerator& rng,
    Deme& deme, Clone& clone, DriverGenotype& driver_genotype, const InputParameters& params) {
    // update event counter for birth
    event_counter.birth++;
    deme.increment(1, params);

    std::vector<int> new_birth_drivers(2,0), new_mig_drivers(2,0);
    std::vector<int> new_mutations = choose_number_mutations(rng, params.mu_driver_birth, params.mu_driver_migration, new_birth_drivers, new_mig_drivers);
    if(!new_mutations[0] && !new_mutations[1]) {
        driver_genotype.increment(1);
        create_clone(deme, driver_genotype, params, event_counter, rng);
        methylation(clone, driver_genotype, params, event_counter, rng);
    }
    else {
        // update event counter for mutations
        event_counter.mutation += new_mutations[0] + new_mutations[1];

        for (int i = 0; i < 2; i++) {
            if(new_birth_drivers[i] || new_mig_drivers[i]) {
                // create clone and perform methylation
                create_clone(deme, driver_genotype, params, event_counter, rng);
                // update driver id
                clones.back().driver_genotype = next_driver_genotype_id;
                create_driver_genotype(clone, driver_genotype);
                // update number of mutations
                driver_genotypes.back().number_of_driver_mutations = new_birth_drivers[i];
                driver_genotypes.back().number_of_migration_mutations = new_mig_drivers[i];
                driver_genotypes.back().set_birth_rate(params, rng);
                driver_genotypes.back().set_migration_rate(params, rng);
            }
            else {
                methylation(clone, driver_genotype, params, event_counter, rng);
            }
        }
        
        if (new_mutations[0] && new_mutations[1]) {
            driver_genotype.increment(-1);
            if (!driver_genotype.immortal) {
                remove_driver_genotype(driver_genotype);
            }
        }
    }
    calculate_deme_birth_rate(deme);
    calculate_deme_migration_rate(deme);
}
// death of cancer cells
void Tumour::cell_death(EventCounter& event_counter, Deme& deme, Clone& clone,
    DriverGenotype& driver_genotype, const InputParameters& params) {
    event_counter.death++;
    deme.increment(-1, params);
    driver_genotype.increment(-1);
    if (!driver_genotype.immortal) {
        remove_driver_genotype(driver_genotype);
    }
    remove_clone(deme, clone);
}

// perform methylation without altering the genotype or clone
void Tumour::methylation(Clone& clone, DriverGenotype& driver_genotype, const InputParameters& params,
    EventCounter& event_counter, RandomNumberGenerator& rng) {
    for (int i = 0; i < clone.meth_array.size(); i++) {
        int condition1 = clone.meth_array[i] == 0 && rng.unif_ran() < params.meth_rate;
        int condition2 = clone.meth_array[i] == 1 && rng.unif_ran() < params.demeth_rate;

        clone.meth_array[i] = clone.meth_array[i] + condition1 - condition2;
        driver_genotype.num_meth += condition1;
        driver_genotype.num_demeth += condition2;
        event_counter.methylation += condition1;
        event_counter.demethylation += condition2;
    }
}

// create new clone upon division and perform methylation - still part of parents' driver genotype
void Tumour::create_clone(Deme& deme, DriverGenotype& parent, const InputParameters& params, 
    EventCounter& event_counter, RandomNumberGenerator& rng) {
    Clone clone(1, deme.identity, next_genotype_id++, parent.driver_identity, deme.clones_list.size());
    methylation(clone, parent, params, event_counter, rng);
    clones.push_back(clone);
    deme.clones_list.push_back(&clone);
}
// create new driver genotype upon division and add to driver genotypes
void Tumour::create_driver_genotype(Clone& clone, DriverGenotype& parent) {
    DriverGenotype driver_genotype(1, clone.driver_genotype, next_driver_genotype_id++, 
        parent.number_of_driver_mutations, parent.number_of_migration_mutations, parent.num_meth,
        parent.num_demeth, 1, parent.birth_rate, parent.migration_rate, gens_elapsed);
    driver_genotypes.push_back(driver_genotype);
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

int Tumour::num_demes() {
    return demes.size();
}

// clone rates
float Tumour::get_clone_birth(Clone& clone) {
    int driver_id = clone.driver_genotype;
    auto it = std::find_if(driver_genotypes.begin(), driver_genotypes.end(), 
                           [driver_id](const DriverGenotype& d) { return d.driver_identity == driver_id; });
    if (it != driver_genotypes.end()) {
        it->birth_rate;
    } else {
        std::cout << "Error: driver genotype not found" << std::endl;
        exit(1);
    }
}

float Tumour::get_clone_migration(Clone& clone) {
    int driver_id = clone.driver_genotype;
    auto it = std::find_if(driver_genotypes.begin(), driver_genotypes.end(), 
                           [driver_id](const DriverGenotype& d) { return d.driver_identity == driver_id; });
    if (it != driver_genotypes.end()) {
        it->migration_rate;
    } else {
        std::cout << "Error: driver genotype not found" << std::endl;
        exit(1);
    }
}

// deme rates
void Tumour::calculate_deme_birth_rate(Deme& deme) {
    float birth_rate = 0;
    for (int i = 0; i < deme.clones_list.size(); i++) {
        int driver_id = deme.clones_list[i]->driver_genotype;
        auto it = std::find_if(driver_genotypes.begin(), driver_genotypes.end(), 
                               [driver_id](const DriverGenotype& g) { return g.driver_identity == driver_id; });
        if (it != driver_genotypes.end()) {
            birth_rate += it->birth_rate;
        } else {
            std::cout << "Error: genotype not found" << std::endl;
            exit(1);
        }
    }
    deme.sum_birth_rates = birth_rate;
}

void Tumour::calculate_deme_migration_rate(Deme& deme) {
    float migration_rate = 0;
    for (int i = 0; i < deme.clones_list.size(); i++) {
        int driver_id = deme.clones_list[i]->driver_genotype;
        auto it = std::find_if(driver_genotypes.begin(), driver_genotypes.end(), 
                               [driver_id](const DriverGenotype& g) { return g.driver_identity == driver_id; });
        if (it != deme.clones_list.end()) {
            migration_rate += it->migration_rate;
        } else {
            std::cout << "Error: genotype not found" << std::endl;
            exit(1);
        }
    }
    deme.sum_migration_rates = migration_rate;
}