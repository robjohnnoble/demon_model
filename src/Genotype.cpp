#include "Objects.hpp"

/////// Genotype
// Constructor
Genotype::Genotype(int population, int parent, int identity, int driverIdentity, int numberOfDriverMutations,
    int numberOfMigrationMutations, int numMeth, int numDemeth, bool immortal, float birthRate, float migrationRate, float originTime)
    : population(population), parent(parent), identity(identity), driver_identity(driverIdentity),
        number_of_driver_mutations(numberOfDriverMutations), number_of_migration_mutations(numberOfMigrationMutations),
        num_meth(numMeth), num_demeth(numDemeth), immortal(immortal), birth_rate(birthRate), migration_rate(migrationRate), origin_time(originTime) {}

// Methods

/////// DriverGenotype
// Constructor
DriverGenotype::DriverGenotype(int population, int parent, int identity, int driverIdentity, int numberOfDriverMutations,
    int numberOfMigrationMutations, int numMeth, int numDemeth, bool immortal, float birthRate, float migrationRate, float originTime)
    : population(population), parent(parent), identity(identity), driver_identity(driverIdentity),
        number_of_driver_mutations(numberOfDriverMutations), number_of_migration_mutations(numberOfMigrationMutations),
        num_meth(numMeth), num_demeth(numDemeth), immortal(immortal), birth_rate(birthRate), migration_rate(migrationRate), origin_time(originTime) {}

// Methods
void DriverGenotype::increment(int increment) {
    population += increment;
    if (population == 0) {
        immortal = false;
    }
}