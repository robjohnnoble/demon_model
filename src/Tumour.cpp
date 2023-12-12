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
    
    // create and shuffle a list of indices
    std::vector<int> indices(parent.population); // vector of indices of clones in parent deme
    std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, 2, ..., parent.population - 1
    std::shuffle(indices.begin(), indices.end(), rng.get_engine()); // mix 'em up real good
    
    // set of indices to move
    std::set<int, std::greater<int> > indices_to_move;
    for (int i = 0; i < num_cells_to_move; ++i) {
        indices_to_move.insert(indices[i]);
    }

    // move the corresponding clones
    for (int i = 0; i < num_cells_to_move; i++) {
        int index = indices[i]; // index of clones in parent deme
        int clone_to_move = clones[parent.clones_list[index]].index; // index of clone in tumour
        std::cout << "clone to move" << clone_to_move << std::endl
            << "clone index in deme " << clones[clone_to_move].index_in_deme << std::endl;
        daughter.clones_list.push_back(clone_to_move);
        // update deme for moved cell
        clones[clone_to_move].deme = daughter.identity;
        // update index in deme for moved cell
        clones[clone_to_move].index_in_deme = daughter.clones_list.size() - 1;
        parent.clones_list.erase(parent.clones_list.begin() + index);
    }
    // update clone indices in parent deme
    int new_population = parent.clones_list.size();
    for (int i = 0; i < new_population; i++) {
        int clone_to_update = clones[parent.clones_list[i]].index;
        clones[clone_to_update].index_in_deme = i;
        std::cout << "clone to update " << clone_to_update << std::endl
            << "clone index in deme " << clones[clone_to_update].index_in_deme << std::endl;
    }
    parent.increment(-num_cells_to_move, params, "deme fission parent");
    daughter.increment(num_cells_to_move, params, "deme fission daughter");
    // std::cout << "successfully moved cells from deme " << parent.identity << " to deme " << daughter.identity << std::endl;
    // for (int i = 0; i < parent.clones_list.size(); i++) {
    //     std::cout << "clone " << i << " in deme " << parent.identity << " has index " << parent.clones_list[i] << std::endl;
    // }
}

// pseudo fission
void Tumour::pseudo_fission(Deme& parent, RandomNumberGenerator& rng, const InputParameters& params) {
    int cells_to_move = rng.stochastic_round(parent.population / 2);
    
    // create and shuffle a list of indices
    int num_clones_in_deme = parent.clones_list.size();
    std::vector<int> indices(num_clones_in_deme);
    std::iota(indices.begin(), indices.end(), 0); // fill with 0, 1, 2, ..., parent.population - 1
    std::shuffle(indices.begin(), indices.end(), rng.get_engine());
    
    // set of indices to remove
    std::set<int, std::greater<int> > indices_to_remove;
    for (int i = 0; i < cells_to_move; i++) {
        indices_to_remove.insert(indices[i]);
    }

    // remove the corresponding clones
    for (int index : indices_to_remove) {
        if (index < parent.clones_list.size()) {
            remove_clone(parent, clones[parent.clones_list[index]]);
        } else {
            // std::cout << "problem when removing clone " << index << "from deme" << parent.identity << std::endl;
        }
    }

    parent.increment(-cells_to_move, params, "pseudo fission");
    calculate_deme_birth_rate(parent);
    calculate_deme_migration_rate(parent);
    parent.set_death_rate(params);
}

// remove clone from deme.clones_list and clones
void Tumour::remove_clone(Deme& deme, Clone& clone) {
    // Handle invalid index, e.g., log an error or throw an exception
    if (clone.index_in_deme >= deme.clones_list.size() || clone.index >= clones.size()) {
        std::cout << "clone index in deme: " << clone.index_in_deme
            << "; clone index: " << clone.index
            << "; chosen deme: " << deme.identity
            << "; population: " << deme.population
            << "; deme clones list size: " << deme.clones_list.size()
            << std::endl;
        return;
    }

    // remove clone from deme
    removeCloneFromDeme(clone.index, deme.identity);
    // remove clone from clones
    if (clone.index < clones.size()) {
        std::swap(clones[clone.index], clones.back());
        clones.pop_back();
        // Update the index of the moved clone
        clones[clone.index].index = clone.index;
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
    // std::cout << "deme population before division " << demes[chosen_deme].clones_list.size() << std::endl;

    std::vector<int> new_birth_drivers(2,0), new_mig_drivers(2,0);
    std::vector<int> new_mutations = choose_number_mutations(rng, params.mu_driver_birth, params.mu_driver_migration, new_birth_drivers, new_mig_drivers);
    if(!new_mutations[0] && !new_mutations[1]) {
        driver_genotypes[clones[chosen_clone].driver_index].increment(1);
        create_clone(clones[chosen_clone], demes[chosen_deme], driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index], params, event_counter, rng);
        demes[chosen_deme].increment(1, params, "cell division");
    }
    else {
        // update event counter for mutations
        event_counter.mutation += new_mutations[0] + new_mutations[1];

        for (int i = 0; i < 2; i++) {
            if(new_birth_drivers[i] || new_mig_drivers[i]) {
                // create clone and perform methylation
                create_clone(clones[chosen_clone], demes[chosen_deme],
                    driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index],
                    params, event_counter, rng);
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
        demes[chosen_deme].increment(1, params, "cell division with mutation");

        if (new_mutations[0] && new_mutations[1]) {
            driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].increment(-1);
            if (!driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].immortal) {
                remove_driver_genotype(driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index]);
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
    driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].increment(-1);

    // if the driver genotype has reached a population of 0, remove it
    if (!driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index].immortal) {
        remove_driver_genotype(driver_genotypes[clones[demes[chosen_deme].clones_list[chosen_clone]].driver_index]);
    }
    // remove clone from deme.clones_list and clones
    remove_clone(demes[chosen_deme], clones[demes[chosen_deme].clones_list[chosen_clone]]);

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
        Clone daughter(1, deme.identity, next_genotype_id++, parent.driver_genotype,
            deme.clones_list.size(), parent.driver_index, clones.size());
        
        // copy methylation array from parent
        daughter.meth_array = parent.meth_array;
        // perform methylation
        methylation(daughter, parent_genotype, params, event_counter, rng);
        // add clone to deme
        deme.clones_list.push_back(daughter.index);
        // add clone to tumour
        clones.push_back(daughter);
}   

// remove clone index from target deme
void Tumour::removeCloneFromDeme(int cloneIndex, int demeIndex) {
    if (demeIndex < demes.size()) {
        Deme& deme = demes[demeIndex];

        // Find the position of the clone's index in Deme::clones_list
        auto it = std::find(deme.clones_list.begin(), deme.clones_list.end(), cloneIndex);

        if (it != deme.clones_list.end()) {
            // Remove the index from Deme::clones_list
            deme.clones_list.erase(it);

            // Update index_in_deme for subsequent clones
            for (size_t i = std::distance(deme.clones_list.begin(), it); i < deme.clones_list.size(); ++i) {
                int updatedCloneIndex = deme.clones_list[i];
                if (updatedCloneIndex < clones.size()) {
                    clones[updatedCloneIndex].index_in_deme = i;
                }
            }
        }
    }
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
    for (int i = 0; i < num_clones_in_deme; i++) {
        int clones_index = deme.clones_list[i];
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
    for (int i = 0; i < num_clones_in_deme; i++) {
        int clones_index = deme.clones_list[i];
        migration_rate += get_clone_migration(clones[clones_index]);
    }
    deme.sum_migration_rates = migration_rate;
}