#include "Objects.hpp"

/////// DriverGenotype
// Constructor
DriverGenotype::DriverGenotype(int population, int parent, int driverIdentity, int numberOfDriverMutations,
    int numberOfMigrationMutations, int numMeth, int numDemeth, bool immortal, float birthRate, float migrationRate, float originTime)
    : population(population), parent(parent), driver_identity(driverIdentity),
        number_of_driver_mutations(numberOfDriverMutations), number_of_migration_mutations(numberOfMigrationMutations),
        num_meth(numMeth), num_demeth(numDemeth), immortal(immortal), birth_rate(birthRate), migration_rate(migrationRate), origin_time(originTime) {}

// Methods
void DriverGenotype::increment(int increment) {
    population += increment;
    if (population == 0) {
        immortal = false;
    }
}

void DriverGenotype::set_birth_rate(const InputParameters& params, RandomNumberGenerator& rng) {
    float advantage = params.s_driver_birth;
    birth_rate = 1;

    if (params.max_relative_birth_rate >= 0) 
        for(int i = 0; i < number_of_driver_mutations; i++) {
            birth_rate = birth_rate * (1 + advantage * (1 - birth_rate / params.max_relative_birth_rate) * rng.exp(1));
        }
    else
        for(int i = 0; i < number_of_driver_mutations; i++) {
            birth_rate = birth_rate * (1 + advantage * rng.exp(1));
        }
    
    if (birth_rate >= params.max_relative_birth_rate + 100) {
        std::cout << "Birth rate exceeds death rate" << std::endl;
        exit(1);
    }
}

void DriverGenotype::set_migration_rate(const InputParameters& params, RandomNumberGenerator& rng) {
    float advantage = params.s_driver_migration;
    migration_rate = params.init_migration_rate;

    if (params.max_relative_migration_rate >= 0) 
        for(int i = 0; i < number_of_migration_mutations; i++) {
            migration_rate = migration_rate * (1 + advantage * (1 - migration_rate / (params.max_relative_migration_rate + params.init_migration_rate)) * rng.exp(1));
        }
    else
        for(int i = 0; i < number_of_migration_mutations; i++) {
            migration_rate = migration_rate * (1 + advantage * rng.exp(1));
        }
}