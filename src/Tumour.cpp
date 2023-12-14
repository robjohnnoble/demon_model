#include "tumour.hpp"

/////// Initialise tumour
void Tumour::initialise(const InputParameters& params,
    const DerivedParameters& d_params) {
    demes.clear();
    genotypes.clear();
    // driver genotypes:
    Genotype firstGenotype(params.init_pop, 0, 0, 0, 0, 0, 0, 0, true, 0);
    genotypes.push_back(firstGenotype);

    // demes:
    Deme firstDeme(d_params.K, "left", 0, 1, 0, params.baseline_death_rate, 1, params.init_migration_rate);
    demes.push_back(firstDeme);
    firstDeme.initialise(params, d_params);

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
            cum_rates[i] += demes[i].sum_rates + cum_rates[i - 1];
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
    const std::set<int>& clones_set = demes[chosen_deme].clones_list;
    std::vector<double> cum_rates;
    double r;

    if (clones_set.size() == 1) {
        return *clones_set.begin();
    } else {
        double rate_sum = 0.0;
        for (int clone_index : clones_set) {
            int driver_id = clones[clone_index].driver_index;
            rate_sum += demes[chosen_deme].death_rate + driver_genotypes[driver_id].birth_rate;
            cum_rates.push_back(rate_sum);
        }
        r = rng.unif_ran() * cum_rates.back();
    }

    if (clones_set.size() == 2) {
        auto it = clones_set.begin();
        return r < cum_rates[0] ? *it : *std::next(it);
    } else {
        auto it = std::lower_bound(cum_rates.begin(), cum_rates.end(), r);
        // Find the corresponding clone index in the set
        auto set_it = std::next(clones_set.begin(), std::distance(cum_rates.begin(), it));
        return *set_it;
    }
}

// choose event type
std::string Tumour::choose_event_type(int chosen_deme, int chosen_clone, RandomNumberGenerator& rng) {
    std::vector<float> cum_rates;
    int ctr = 0;
    float res;

    cum_rates.push_back(get_clone_birth(clones[chosen_clone]));
    ctr++;

    cum_rates.push_back(cum_rates[ctr - 1] + demes[chosen_deme].death_rate);
    ctr++;

    cum_rates.push_back(cum_rates[ctr - 1] + demes[chosen_deme].sum_migration_rates);

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

// indicator for fission at birth event
bool Tumour::fission_ready(int chosen_deme, RandomNumberGenerator& rng, bool birth) {
    if (!birth) return demes[chosen_deme].population >= demes[chosen_deme].K;
    else {
        return demes[chosen_deme].population >= demes[chosen_deme].K && 
            demes[chosen_deme].sum_migration_rates >= rng.unif_ran();
    }
}

// perform deme fission
void Tumour::deme_fission(int chosen_deme, EventCounter& event_counter, RandomNumberGenerator& rng, const InputParameters& params){
    int new_deme_id = static_cast<int>(demes.size());
    if (next_fission == fission_times[0]) {
        // first fission - sets up left and right sides of the tumour
        event_counter.fission++;
        // create new deme
        Deme new_deme(demes[chosen_deme].K, std::string("right"), new_deme_id, 0, demes[chosen_deme].fissions + 1, 0, 0, 0);
        new_deme.clones_list.clear();
        // move cells into new deme
        move_cells(demes[chosen_deme], new_deme, rng, params);
        // add new deme to tumour
        demes.push_back(new_deme);
        // update next fission time
        next_fission = fission_times[1];
    } else if (gens_elapsed >= next_fission && new_deme_id < 8) {
        // subsequent fissions
        event_counter.fission++;
        Deme new_deme(demes[chosen_deme].K, demes[chosen_deme].side, new_deme_id, 0, demes[chosen_deme].fissions + 1, 0, 0, 0);
        new_deme.clones_list.clear();
        move_cells(demes[chosen_deme], new_deme, rng, params);
        demes.push_back(new_deme);
        // update next fission time
        if(next_fission != fission_times[6]) {
            int tmp = demes.size() - 1;
            next_fission = fission_times[tmp];
        }
        // if the last fission happened, set next fission to be after max_generations
        else next_fission = params.max_generations + 1;
    } else {
        // pseudo fission - kills half the population in a deme without creating a new deme
        event_counter.fission++;
        // counts as a fission for the deme
        demes[chosen_deme].fissions++;
        pseudo_fission(demes[chosen_deme], rng, params);
    }
}
// move cells into new deme after deme fission
void Tumour::move_cells(Deme& parent, Deme& daughter, RandomNumberGenerator& rng, const InputParameters& params) {
    int num_cells_to_move = rng.stochastic_round(parent.population / 2);

    // Create a vector from the set for easy random access
    std::vector<int> parent_indices(parent.clones_list.begin(), parent.clones_list.end());

    // Shuffle the vector to randomize which cells are moved
    std::shuffle(parent_indices.begin(), parent_indices.end(), rng.get_engine());

    // Move the specified number of cells
    for (int i = 0; i < num_cells_to_move && i < parent_indices.size(); ++i) {
        int clone_index = parent_indices[i];

        // Move clone from parent to daughter
        parent.clones_list.erase(clone_index);
        daughter.clones_list.insert(clone_index);

        // Update the deme field of the moved clone
        clones[clone_index].deme = daughter.identity;
    }

    // Update the populations of parent and daughter demes
    parent.increment(-num_cells_to_move, params, "deme fission");
    daughter.increment(num_cells_to_move, params, "deme fission");
}

// pseudo fission
void Tumour::pseudo_fission(Deme& parent, RandomNumberGenerator& rng, const InputParameters& params) {
    int cells_to_remove = rng.stochastic_round(parent.population / 2);
    
    // Create a vector from the set for easy random access
    std::vector<int> parent_indices(parent.clones_list.begin(), parent.clones_list.end());

    // Shuffle the vector to randomize which cells are moved
    std::shuffle(parent_indices.begin(), parent_indices.end(), rng.get_engine());

    // Remove the specified number of cells
    for (int i = 0; i < cells_to_remove && i < parent_indices.size(); ++i) {
        int clone_index = parent_indices[i];
        if (parent.clones_list.find(clone_index) == parent.clones_list.end()) {
            std::cout << "Clone not found in parent" << std::endl;
            exit(1);
        }
        // Remove clone from parent
        parent.clones_list.erase(clone_index);
        // Remove clone from tumour
        remove_clone(clones[clone_index]);
    }
    parent.increment(-cells_to_remove, params, "pseudo fission");
}

// remove clone from tumour.clones
void Tumour::remove_clone(Clone& clone) {
    std::swap(clones[clone.index], clones.back());
    clones.pop_back();
    // update all demes' clones lists
    for (int i = 0; i < demes.size(); i++) {
        demes[i].clones_list.clear();
    }
    for (int i = 0; i < clones.size(); i++) {
        clones[i].index = i;
        int deme_index = clones[i].deme;
        demes[deme_index].clones_list.insert(i);
    }
}
// remove clone index from target deme
// void Tumour::remove_clone_from_deme(int cloneIndex, int demeIndex) {
//     demes[demeIndex].clones_list.erase(cloneIndex);
// }

// remove driver genotype from driver genotypes
void Tumour::remove_driver_genotype(DriverGenotype& driver_genotype) {
    int old_index = driver_genotype.index;
    std::swap(driver_genotypes[old_index], driver_genotypes.back());
    driver_genotypes.pop_back();
    for (int i = 0; i < driver_genotypes.size(); i++) {
        driver_genotypes[i].index = i;
    }
    for (int i = 0; i < clones.size(); i++) {
        if (clones[i].driver_index == driver_genotypes.size()) {
            clones[i].driver_index = old_index;
        }
    }
}


// division of cancer cells
void Tumour::cell_division(EventCounter& event_counter, RandomNumberGenerator& rng,
    int chosen_deme, int chosen_clone, const InputParameters& params) {
    // update event counter for birth
    event_counter.birth++;
    // std::cout << "deme population before division " << demes[chosen_deme].clones_list.size() << std::endl;

    std::vector<int> new_birth_drivers(2,0), new_mig_drivers(2,0);
    std::vector<int> new_mutations = choose_number_mutations(rng, params.mu_driver_birth, params.mu_driver_migration, new_birth_drivers, new_mig_drivers);
    if(!new_mutations[0] && !new_mutations[1]) {
        driver_genotypes[clones[chosen_clone].driver_index].increment(1);
        create_clone(clones[chosen_clone], demes[chosen_deme], driver_genotypes[clones[chosen_clone].driver_index], params, event_counter, rng);
        demes[chosen_deme].increment(1, params, "cell division");
    }
    else {
        // update event counter for mutations
        event_counter.mutation += new_mutations[0] + new_mutations[1];

        for (int i = 0; i < 2; i++) {
            if(new_birth_drivers[i] || new_mig_drivers[i]) {
                // create clone and perform methylation
                create_clone(clones[chosen_clone], demes[chosen_deme],
                    driver_genotypes[clones[chosen_clone].driver_index],
                    params, event_counter, rng);
                // update driver id
                clones.back().driver_genotype = next_driver_genotype_id;
                create_driver_genotype(clones[chosen_clone], driver_genotypes[clones[chosen_clone].driver_index]);
                clones.back().driver_index = driver_genotypes.size() - 1;
                // update number of mutations
                driver_genotypes.back().number_of_driver_mutations += new_birth_drivers[i];
                driver_genotypes.back().number_of_migration_mutations += new_mig_drivers[i];
                driver_genotypes.back().set_birth_rate(params, rng);
                driver_genotypes.back().set_migration_rate(params, rng);
            }
            else {
                methylation(clones[chosen_clone], driver_genotypes[clones[chosen_clone].driver_index], params, event_counter, rng);
            }
        }
        demes[chosen_deme].increment(1, params, "cell division with mutation");

        if (new_mutations[0] && new_mutations[1]) {
            driver_genotypes[clones[chosen_clone].driver_index].increment(-1);
            if (!driver_genotypes[clones[chosen_clone].driver_index].immortal) {
                remove_driver_genotype(driver_genotypes[clones[chosen_clone].driver_index]);
            }
        }
    }
    calculate_deme_birth_rate(demes[chosen_deme]);
    calculate_deme_migration_rate(demes[chosen_deme]);
    // std::cout << "deme" << demes[chosen_deme].identity << " population after division "
    //     << demes[chosen_deme].clones_list.size() << std::endl;
}
// death of cancer cells
void Tumour::cell_death(EventCounter& event_counter, int chosen_deme, int chosen_clone, const InputParameters& params) {
    // update event counter for death
    event_counter.death++;

    // update driver genotype population
    driver_genotypes[clones[chosen_clone].driver_index].increment(-1);

    // if the driver genotype has reached a population of 0, remove it
    if (!driver_genotypes[clones[chosen_clone].driver_index].immortal) {
        remove_driver_genotype(driver_genotypes[clones[chosen_clone].driver_index]);
    }
    // remove clone from deme.clones_list and clones
    remove_clone(clones[chosen_clone]);
    demes[chosen_deme].increment(-1, params, "cell death");
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

void Tumour::calculate_average_array(int deme_index, const DerivedParameters& d_params) {
    std::vector<float> res(d_params.fcpgs, 0);
    demes[deme_index].avg_meth_array = std::vector<float>(d_params.fcpgs, 0);

    if (!demes[deme_index].clones_list.size()) std::cout << "No clones in deme" << std::endl;
    for (int clone_index : demes[deme_index].clones_list) {
        for (int j = 0; j < d_params.fcpgs; j++) {
            res[j] += clones[clone_index].meth_array[j];
        }
    }

    for (int i = 0; i < d_params.fcpgs; i++) {
        demes[deme_index].avg_meth_array[i] = res[i] / static_cast<float>(demes[deme_index].population);
    }
}

// create new clone upon division and perform methylation - still part of parents' driver genotype
void Tumour::create_clone(const Clone& parent, Deme& deme, DriverGenotype& parent_genotype, const InputParameters& params, 
    EventCounter& event_counter, RandomNumberGenerator& rng) {
        int next_clone_index = clones.size();
        Clone daughter(1, deme.identity, next_genotype_id++, parent.driver_genotype,
            parent.driver_index, next_clone_index);
        
        // copy methylation array from parent and perform methylation
        daughter.meth_array = parent.meth_array;
        methylation(daughter, parent_genotype, params, event_counter, rng);

        // add clone to deme
        deme.clones_list.insert(daughter.index);
        // add clone to tumour
        clones.push_back(daughter);
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
    sum_death_rates = 0;
    sum_birth_rates = 0;
    sum_migration_rates = 0;
    for (int i = 0; i < demes.size(); i++) {
        demes[i].calculate_sum_of_rates();
        sum_death_rates += demes[i].death_rate * demes[i].population;
        sum_birth_rates += demes[i].sum_birth_rates;
        sum_migration_rates += demes[i].sum_migration_rates;
    }
}

double Tumour::sum_of_all_rates() {
    calculate_sums_of_rates();
    return sum_birth_rates + sum_death_rates + sum_migration_rates;
}

// update time
void Tumour::update_time(RandomNumberGenerator& rng, const InputParameters& params) {
    float temp_sum = sum_of_all_rates();
    // std::cout << "Sum of all rates: " << temp_sum << std::endl;

    float gens_added = rng.exp(1) / temp_sum;
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
    int num_clones_in_deme = deme.clones_list.size();
    if(!num_clones_in_deme) {
        std::cout << "No clones in deme" << std::endl;
        return;
    }
    for (int clones_index : deme.clones_list) {
        birth_rate += get_clone_birth(clones[clones_index]);
    }
    deme.sum_birth_rates = birth_rate;
}

void Tumour::calculate_deme_migration_rate(Deme& deme) {
    float migration_rate = 0;
    int num_clones_in_deme = deme.clones_list.size();
    if(!num_clones_in_deme) {
        std::cout << "No clones in deme" << std::endl;
        return;
    }
    for (int clones_index : deme.clones_list) {
        migration_rate += get_clone_migration(clones[clones_index]);
    }
    deme.sum_migration_rates = migration_rate;
}

// update all rates
void Tumour::calculate_all_rates(const InputParameters& params, const DerivedParameters& d_params) {
    for (int i = 0; i < demes.size(); i++) {
        calculate_deme_birth_rate(demes[i]);
        calculate_deme_migration_rate(demes[i]);
        demes[i].set_death_rate(params);
        demes[i].calculate_sum_of_rates();
    }
}