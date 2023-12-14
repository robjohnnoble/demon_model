#ifndef TUMOUR_HPP
#define TUMOUR_HPP

#include "parameters.hpp"
#include "distributions.hpp"
#include "deme.hpp"
#include "genotype.hpp"
#include "output.hpp"
#include "macros.hpp"
#include <set>
#include <functional>
#include <vector>
#include <iostream>

class Tumour {
    private:
    // cell containers
    std::vector<Deme> demes;
    std::vector<Genotype> genotypes;
    // cell tracking
    int next_cell_id;
    // aux variables
    int next_genotype_id = 1;
    int next_driver_genotype_id = 1;

    // temporal variables
    float gens_elapsed = 0;
    float output_timer = 0;
    std::vector<int> fission_times;
    int next_fission = 0;    

    // sums of rates
    double sum_birth_rates = 0;
    double sum_death_rates = 0;
    double sum_migration_rates = 0;

    public:
    // initialise the tumour from input parameters
    void initialise(const InputParameters& params, const DerivedParameters& d_params);
    // sum of different rates in the tumour
    void calculate_sums_of_rates();
    double sum_of_all_rates();

    // deme fission
    bool fission_ready(int chosen_deme, RandomNumberGenerator& rng, bool birth);
    void deme_fission(int chosen_deme, EventCounter& event_counter, RandomNumberGenerator& rng, const InputParameters& params);
    void remove_cell(Cell& cell);
    void remove_cell_from_deme(int cellIndex, int demeIndex);
    void move_cells(Deme& parent, Deme& daughter, RandomNumberGenerator& rng, const InputParameters& params);
    void pseudo_fission(Deme& deme, RandomNumberGenerator& rng, const InputParameters& params);

    // cell division
    void cell_division(EventCounter& event_counter, RandomNumberGenerator& rng,
        int chosen_deme, int chosen_cell, const InputParameters& params);
    // cell death
    void cell_death(EventCounter& event_counter, int chosen_deme, int chosen_cell, const InputParameters& params);
    void remove_driver_genotype(DriverGenotype& driver_genotype);
    // create new genotype upon division
    void create_cell(const Cell& parent, Deme& deme, DriverGenotype& parent_genotype, const InputParameters& params,
        EventCounter& event_counter, RandomNumberGenerator& rng);
    // create new driver genotype upon division and add to driver genotypes
    void create_driver_genotype(const Cell& cell, DriverGenotype& parent);
    // choose number of mutations upon division
    std::vector<int> choose_number_mutations(RandomNumberGenerator& rng, float mu_driver_birth, float mu_driver_migration, 
        std::vector<int>& new_birth_drivers, std::vector<int>& new_mig_drivers);
    // methylation
    void methylation(Cell& cell, DriverGenotype& driver_genotype, const InputParameters& params,
        EventCounter& event_counter, RandomNumberGenerator& rng);
    void calculate_average_array(int deme_index, const DerivedParameters& d_params);

    // cell rates
    float get_cell_birth(const Cell& cell);
    float get_cell_migration(const Cell& cell);
    
    // deme rates
    void calculate_deme_birth_rate(Deme& deme);
    void calculate_deme_migration_rate(Deme& deme);
    void calculate_all_rates(const InputParameters& params, const DerivedParameters& d_params);
    
    // choose deme, cell and event type
    int choose_deme(RandomNumberGenerator& rng);
    int choose_cell(int chosen_deme, RandomNumberGenerator& rng);
    std::string choose_event_type(int chosen_deme, int chosen_cell, RandomNumberGenerator& rng);

    // numbers of relevant things
    int num_cells();
    int num_cells();
    int num_driver_genotypes();
    int num_demes();
    
    // update time
    void update_time(RandomNumberGenerator& rng, const InputParameters& params);
    // check time
    float check_time();
    
    // output
    void print_to_screen(const InputParameters& params, const DerivedParameters& d_params, const EventCounter& event_counter);
    //void final_output();
};

#endif // TUMOUR_HPP