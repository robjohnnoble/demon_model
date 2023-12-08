#ifndef OBJECTS_HPP
#define OBJECTS_HPP

#include "Parameters.hpp"
#include "Distributions.hpp"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

class Deme {
    public:
    // Constructor
    Deme(int K, std::string side, int identity, int population, int fissions, float deathRate, float migrationModifier, float sumBirthRates, float sumMigrationRates);

    // Methods
    void calculate_sum_of_rates();
    void calculate_average_array(const DerivedParameters& d_params, std::vector<Clone>& clones);
    void increment(int increment, const InputParameters& params);
    void set_death_rate(const InputParameters& params);
    void remove_clone(Clone& clone);

    // Properties
    int K; // carrying capacity of the deme

    std::string side; // Left or right

    int identity; // Identity of the deme
    int population; // Number of cancer cells in the deme
    std::vector<Clone*> clones_list; // List of clones in the deme
    int fissions; // fissions since the initial deme

    float death_rate; // Death rate of the deme
    float migration_modifier; // factor by which migration rate is multiplied (depends on what cells are in the deme)
    std::vector<float> avg_meth_array; // Average methylation array of the deme
    float sum_birth_rates; // Sum of birth rates of the clones in the deme
    float sum_migration_rates; // Sum of migration rates of the clones in the deme
    float sum_rates; // sum of all rates in the deme
};

class DriverGenotype {
    public:
    // Properties
    int population;
    int parent; // parent's unique ID
    int driver_identity; // unique ID of the driver genotype
    int number_of_driver_mutations; // number of driver mutations
    int number_of_migration_mutations; // number of migration mutations
    int num_meth; // number of methylation events since initial array
    int num_demeth; // number of demethylation events since initial array

    bool immortal; // whether genotype record can be overwritten

    float birth_rate; // birth rate of the genotype
    float migration_rate; // migration rate of the genotype
    float origin_time; // generation in which the genotype was created

    int index;

    // Constructor
    DriverGenotype(int population, int parent, int driverIdentity, int numberOfDriverMutations, int numberOfMigrationMutations, int numMeth, int numDemeth, bool immortal, float birthRate, float migrationRate, float originTime);

    // Methods
    void increment(int increment);
    void set_birth_rate(const InputParameters& params, RandomNumberGenerator& rng);
    void set_migration_rate(const InputParameters& params, RandomNumberGenerator& rng);
};

class Clone {
    public:
    // Properties
    int population; // Number of cancer cells
    int deme; // Index of the deme in which the clone is located
    int genotype; // Identity of the clone
    int driver_genotype; // Driver genotype of the clone
    int driver_index; // Index of the driver genotype in the tumour
    // TO DO: add driver index handling in relevant functions

    int index_in_deme; // Index of the clone within its deme
    int index; // index in tumour

    std::vector<int> meth_array; // fCpG array of the clone

    // Constructor
    Clone(int population, int deme, int genotype, int driverGenotype, int indexInDeme);

    // Methods
    void initial_array(const InputParameters& params, const DerivedParameters& d_params, RandomNumberGenerator& rng);
};

#endif // OBJECTS_HPP