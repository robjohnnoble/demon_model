#include "Tumour.hpp"

// initialise tumour from input parameters
void Tumour::initialise(const InputParameters& params,
    const DerivedParameters& d_params, RandomNumberGenerator& rng) {
    // driver genotypes:
    DriverGenotype driver_genotype(params.init_pop, 0, 0, 0, 0, 0, 0, true, 1, params.init_migration_rate, 0);
    driver_genotype.index = 0;
    driver_genotypes.push_back(driver_genotype);

    // clones:
    Clone clone(1, 0, 0, 0, 0, 0, 0);
    clone.initial_array(params, d_params, rng);
    clones.push_back(clone);
    // for (int i = 0; i < d_params.fcpgs; i++)
    //     std::cout << clone.meth_array[i];
    // std::cout << std::endl;

    // demes:
    Deme deme(d_params.K, "left", 0, 1, 0, params.baseline_death_rate, 1, params.init_migration_rate);
    deme.clones_list.push_back(0);
    deme.calculate_sum_of_rates();
    calculate_average_array(deme, d_params);
    demes.push_back(deme);

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
    int res;

    if (cum_rates.size() == 1) {
        res = 0;
        calculate_deme_birth_rate(demes[res]);
        calculate_deme_migration_rate(demes[res]);
        return res;
    } else {
        cum_rates[0] = demes[0].sum_rates;
        for (int i = 1; i < cum_rates.size(); i++) {
            cum_rates[i] = demes[i].sum_rates + cum_rates[i - 1];
        }
        r = rng.unif_ran() * cum_rates.back();
    }
    if (cum_rates.size() == 2) {
        res = r < cum_rates[0] ? 0 : 1;
        calculate_deme_birth_rate(demes[res]);
        calculate_deme_migration_rate(demes[res]);
        return res;
    } else {
        res = std::lower_bound(cum_rates.begin(), cum_rates.end(), r) - cum_rates.begin();
        calculate_deme_birth_rate(demes[res]);
        calculate_deme_migration_rate(demes[res]);
        return res;
    }
}
// choose clone
int Tumour::choose_clone(int chosen_deme, RandomNumberGenerator& rng) {
    std::vector<double> cum_rates(demes[chosen_deme].clones_list.size());
    double r;

    if (demes[chosen_deme].clones_list.size() == 1) {
        return 0;
    } else {
        int driver_id = clones[demes[chosen_deme].clones_list[0]].driver_index;
        cum_rates[0] = demes[chosen_deme].death_rate + driver_genotypes[driver_id].birth_rate;
        for (int i = 1; i < demes[chosen_deme].clones_list.size(); i++) {
            driver_id = clones[demes[chosen_deme].clones_list[i]].driver_index;
            cum_rates[i] = cum_rates[i-1] + demes[chosen_deme].death_rate + driver_genotypes[driver_id].birth_rate;
        }
        r = rng.unif_ran() * cum_rates.back();
    }

    if (demes[chosen_deme].clones_list.size() == 2) {
        return r < cum_rates[0] ? 0 : 1;
    } else {
        return std::lower_bound(cum_rates.begin(), cum_rates.end(), r) - cum_rates.begin();
    }
}
// choose event type
std::string Tumour::choose_event_type(int chosen_deme, int chosen_clone, RandomNumberGenerator& rng) {
    std::vector<float> cum_rates;
    int ctr = 0;
    float res;

    cum_rates.push_back(get_clone_birth(clones[demes[chosen_deme].clones_list[chosen_clone]]));
    ctr++;

    cum_rates.push_back(cum_rates[ctr - 1] + demes[chosen_deme].death_rate);
    ctr++;

    cum_rates.push_back(cum_rates[ctr - 1] + get_clone_migration(clones[demes[chosen_deme].clones_list[chosen_clone]]) * demes[chosen_deme].sum_migration_rates);

    res = rng.unif_ran() * cum_rates.back();
    if(res < cum_rates[0]) {
        return "birth";
    }
    else if(res < cum_rates[1]) {
        return "death";
    }
    else {
        return "fission";
    }
}

// deme fission
bool Tumour::fission_ready(int chosen_deme, RandomNumberGenerator& rng, bool birth) {
    if (!birth) return demes[chosen_deme].population >= demes[chosen_deme].K;
    else {
        return demes[chosen_deme].population >= demes[chosen_deme].K && 
            demes[chosen_deme].sum_migration_rates >= rng.unif_ran();
    }
}
void Tumour::deme_fission(int chosen_deme, EventCounter& event_counter, RandomNumberGenerator& rng, const InputParameters& params){
    int new_deme_id = static_cast<int>(demes.size());
    if (next_fission == fission_times[0]) {
        // first fission - sets up left and right sides of the tumour
        event_counter.fission++;
        Deme new_deme(demes[chosen_deme].K, std::string("right"), new_deme_id, 0, demes[chosen_deme].fissions + 1, 0, 0, 0);
        move_cells(demes[chosen_deme], new_deme, rng, params);
        demes.push_back(new_deme);
        next_fission = fission_times[1];
    } else if (gens_elapsed >= next_fission) {
        // subsequent fissions
        event_counter.fission++;
        Deme new_deme(demes[chosen_deme].K, demes[chosen_deme].side, new_deme_id, 0, demes[chosen_deme].fissions + 1, 0, 0, 0);
        move_cells(demes[chosen_deme], new_deme, rng, params);
        demes.push_back(new_deme);
        if(next_fission != fission_times.back()) next_fission = fission_times[ - 1];
        else next_fission = params.max_generations + 1;
    } else {
        // pseudo fission - kills half the population in a deme without creating a new deme
        event_counter.fission++;
        demes[chosen_deme].fissions++;
        pseudo_fission(demes[chosen_deme], rng, params);
    }
}
// move cells after deme fission
void Tumour::move_cells(Deme& parent, Deme& daughter, RandomNumberGenerator& rng, const InputParameters& params) {
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
        daughter.clones_list.push_back(parent.clones_list[index]);
        // update deme for moved cell
        clones[daughter.clones_list.back()].deme = daughter.identity;
        // update index in deme for moved cell
        clones[daughter.clones_list.back()].index_in_deme = daughter.clones_list.size() - 1;
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
        remove_clone(parent, clones[parent.clones_list[index]]);
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
    }
}


// division of cancer cells
void Tumour::cell_division(EventCounter& event_counter, RandomNumberGenerator& rng,
    int chosen_deme, int chosen_clone, const InputParameters& params) {
    // update event counter for birth
    event_counter.birth++;
    demes[chosen_deme].increment(1, params);

    std::vector<int> new_birth_drivers(2,0), new_mig_drivers(2,0);
    std::vector<int> new_mutations = choose_number_mutations(rng, params.mu_driver_birth, params.mu_driver_migration, new_birth_drivers, new_mig_drivers);
    if(!new_mutations[0] && !new_mutations[1]) {
        driver_genotypes[clones[chosen_clone].driver_index].increment(1);
        create_clone(clones[chosen_clone], demes[chosen_deme], driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index], params, event_counter, rng);
    }
    else {
        // update event counter for mutations
        event_counter.mutation += new_mutations[0] + new_mutations[1];

        for (int i = 0; i < 2; i++) {
            if(new_birth_drivers[i] || new_mig_drivers[i]) {
                // create clone and perform methylation
                create_clone(clones[chosen_clone], demes[chosen_deme], driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index], params, event_counter, rng);
                // update driver id
                clones.back().driver_genotype = next_driver_genotype_id;
                create_driver_genotype(clones[demes[chosen_deme].clones_list[chosen_clone]], driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index]);
                // update number of mutations
                driver_genotypes.back().number_of_driver_mutations = new_birth_drivers[i];
                driver_genotypes.back().number_of_migration_mutations = new_mig_drivers[i];
                driver_genotypes.back().set_birth_rate(params, rng);
                driver_genotypes.back().set_migration_rate(params, rng);
            }
            else {
                methylation(clones[demes[chosen_deme].clones_list[chosen_clone]], driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index], params, event_counter, rng);
            }
        }
        
        if (new_mutations[0] && new_mutations[1]) {
            driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].increment(-1);
            if (!driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].immortal) {
                remove_driver_genotype(driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index]);
            }
        }
    }
    calculate_deme_birth_rate(demes[chosen_deme]);
    calculate_deme_migration_rate(demes[chosen_deme]);
}
// death of cancer cells
void Tumour::cell_death(EventCounter& event_counter, int chosen_deme, int chosen_clone, const InputParameters& params) {
    event_counter.death++;
    demes[chosen_deme].increment(-1, params);
    driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].increment(-1);
    if (!driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].immortal) {
        remove_driver_genotype(driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index]);
    }
    remove_clone(demes[chosen_deme], clones[demes[chosen_deme].clones_list[chosen_clone]]);
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

void Tumour::calculate_average_array(Deme& deme, const DerivedParameters& d_params) {
    std::vector<float> res(d_params.fcpgs, 0);
    deme.avg_meth_array = std::vector<float>(d_params.fcpgs, 0);

    if (!deme.clones_list.size()) std::cout << "No clones in deme" << std::endl;
    for (int i = 0; i < deme.clones_list.size(); i++) {
        for (int j = 0; j < d_params.fcpgs; j++) {
            res[j] += clones[deme.clones_list[i]].meth_array[j];// * clones_list[i]->population;
        }
    }

    for (int i = 0; i < d_params.fcpgs; i++) {
        deme.avg_meth_array[i] = res[i] / static_cast<float>(deme.population);
    }
}

// create new clone upon division and perform methylation - still part of parents' driver genotype
void Tumour::create_clone(const Clone& parent, Deme& deme, DriverGenotype& parent_genotype, const InputParameters& params, 
    EventCounter& event_counter, RandomNumberGenerator& rng) {
    Clone clone(1, deme.identity, next_genotype_id++, parent_genotype.driver_identity, deme.clones_list.size(), parent.driver_index, clones.size());
    clone.driver_index = parent.driver_index;
    // copy meth array into new clone
    clone.meth_array.clear();
    for (int i = 0; i < parent.meth_array.size(); i++) {
        clone.meth_array.push_back(parent.meth_array[i]);
    }
    // perform methylation
    methylation(clone, parent_genotype, params, event_counter, rng);
    clones.push_back(clone);
    deme.clones_list.push_back(clones.back().index);
}
// create new driver genotype upon division and add to driver genotypes
void Tumour::create_driver_genotype(const Clone& clone, DriverGenotype& parent) {
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
        res += demes[i].clones_list.size();
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
float Tumour::get_clone_birth(const Clone& clone) {
    return driver_genotypes[clone.driver_index].birth_rate;
}

float Tumour::get_clone_migration(const Clone& clone) {
    return driver_genotypes[clone.driver_index].migration_rate;
}

// deme rates
void Tumour::calculate_deme_birth_rate(Deme& deme) {
    float birth_rate = 0;
    for (int i = 0; i < deme.clones_list.size(); i++) {
        birth_rate += get_clone_birth(clones[deme.clones_list[i]]);
    }
    deme.sum_birth_rates = birth_rate;
}

void Tumour::calculate_deme_migration_rate(Deme& deme) {
    float migration_rate = 0;
    for (int i = 0; i < deme.clones_list.size(); i++) {
        migration_rate += get_clone_migration(clones[deme.clones_list[i]]);
    }
    deme.sum_migration_rates = migration_rate;
}