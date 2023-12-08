#ifndef TUMOUR_HPP
#define TUMOUR_HPP

#include "Parameters.hpp"
#include "Distributions.hpp"
#include "Objects.hpp"
#include "Macros.hpp"
#include <vector>
#include <iostream>

class Tumour {
    private:
    // objects
    std::vector<Deme> demes;
    std::vector<Clone> clones;
    std::vector<DriverGenotype> driver_genotypes;
    // aux variables
    int next_genotype_id = 1;
    int next_driver_genotype_id = 1;

    // temporal variables
    float gens_elapsed;
    float output_timer;
    std::vector<int> fission_times;
    int next_fission;    

    // sums of rates
    double sum_birth_rates;
    double sum_death_rates;
    double sum_migration_rates;

    public:
    // initialise the tumour from input parameters
    void initialise(const InputParameters& params, const DerivedParameters& d_params, RandomNumberGenerator& rng);
    // sum of different rates in the tumour
    void calculate_sums_of_rates();
    double sum_of_all_rates();

    // deme fission
    void deme_fission(Deme& deme, EventCounter& event_counter, RandomNumberGenerator& rng, const InputParameters& params);
    void remove_clone(Deme& deme, Clone& clone);
    void move_cells(Deme& parent, Deme& daughter, RandomNumberGenerator& rng, const InputParameters& params);
    void pseudo_fission(Deme& deme, RandomNumberGenerator& rng, const InputParameters& params);

    // cell division
    void cell_division(EventCounter& event_counter, RandomNumberGenerator& rng,
        Deme& deme, Clone& clone, DriverGenotype& driver_genotype, const InputParameters& params);
    // cell death
    void Tumour::cell_death(EventCounter& event_counter, Deme& deme, Clone& clone,
    DriverGenotype& driver_genotype, const InputParameters& params);
    void Tumour::remove_driver_genotype(DriverGenotype& driver_genotype);
    // create new genotype upon division
    void create_clone(Deme& deme, DriverGenotype& parent, const InputParameters& params,
        EventCounter& event_counter, RandomNumberGenerator& rng);
    // create new driver genotype upon division and add to driver genotypes
    void create_driver_genotype(Clone& clone, DriverGenotype& parent);
    // choose number of mutations upon division
    std::vector<int> Tumour::choose_number_mutations(RandomNumberGenerator& rng, float mu_driver_birth, float mu_driver_migration, 
        std::vector<int>& new_birth_drivers, std::vector<int>& new_mig_drivers);
    // methylation
    void methylation(Clone& clone, DriverGenotype& driver_genotype, const InputParameters& params,
        EventCounter& event_counter, RandomNumberGenerator& rng);

    // clone rates
    float get_clone_birth(Clone& clone);
    float get_clone_migration(Clone& clone);
    
    // deme rates
    void calculate_deme_birth_rate(Deme& deme);
    void calculate_deme_migration_rate(Deme& deme);
    
    // choose deme, clone and event type
    int Tumour::choose_deme(RandomNumberGenerator& rng);
    int Tumour::choose_clone(int chosen_deme, RandomNumberGenerator& rng);
    std::string Tumour::choose_event_type();

    // numbers of relevant things
    int num_cells();
    int num_clones();
    int num_driver_genotypes();
    int num_demes();
    
    // update time
    void update_time(RandomNumberGenerator& rng, const InputParameters& params);
    // check time
    float check_time();
    
    // output
    void print_to_screen(Tumour& tumour, const InputParameters& params, const DerivedParameters& d_params, const EventCounter& event_counter);
    void final_output();
};

#endif // TUMOUR_HPP