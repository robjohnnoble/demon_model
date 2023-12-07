#ifndef TUMOUR_HPP
#define TUMOUR_HPP

#include "Parameters.hpp"
#include "Distributions.hpp"
#include "Objects.hpp"
#include <vector>

class Tumour {
    private:
    // objects
    std::vector<Deme> demes;
    std::vector<Clone> clones;
    std::vector<Genotype> genotypes;
    std::vector<DriverGenotype> driver_genotypes;
    // aux variables
    int next_genotype_id = 1;

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
    void deme_fission();

    // cell division
    void cell_division(int parent_clone, EventCounter& event_counter, RandomNumberGenerator& rng,
        Deme& deme, Clone& clone, Genotype& genotype, DriverGenotype& driver_genotype, const InputParameters& params);
    // create new genotype upon division
    void create_genotype(Genotype& parent);
    // choose number of mutations upon division
    std::vector<int> Tumour::choose_number_mutations(RandomNumberGenerator& rng, float mu_driver_birth, float mu_driver_migration, 
        std::vector<int>& new_birth_drivers, std::vector<int>& new_mig_drivers);

    // numbers of relevant things
    int num_cells();
    int num_clones();
    int num_driver_genotypes();
    int num_genotypes();
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