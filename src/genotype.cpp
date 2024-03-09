#include "genotype.hpp"

/////// Constructor
Genotype::Genotype(int parent, int identity, int numBirthMut, int numMigMut,
                   float birthRate, float migrationRate, float originTime,
                   const InputParameters &params)
    : parent(parent), identity(identity), numBirthMut(numBirthMut),
      numMigMut(numMigMut), birthRate(birthRate), migrationRate(migrationRate),
      muDriverBirth(params.mu_driver_birth),
      muDriverMig(params.mu_driver_migration),
      inputMigRate(params.init_migration_rate),
      maxRelBirth(params.max_relative_birth_rate),
      maxRelMig(params.max_relative_migration_rate),
      sDriverBirth(params.s_driver_birth),
      sDriverMig(params.s_driver_migration), originTime(originTime) {}

/////// Setters
// birth rate
void Genotype::setBirthRate() {
    birthRate = 1;

    if (maxRelBirth >= 0)
        for(int i = 0; i < numBirthMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            birthRate = birthRate * (1 + sDriverBirth * (1 - birthRate / maxRelBirth) * rnd);
        }
    else
        for(int i = 0; i < numBirthMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            birthRate = birthRate * (1 + sDriverBirth * rnd);
        }

    if (birthRate >= maxRelBirth + 100) {
        std::cout << "ERROR: Birth rate at carrying capacity exceeds death rate." << std::endl;
        exit(1);
    }
}
// set migration rate
void Genotype::setMigrationRate() {
    migrationRate = inputMigRate;

    if (maxRelMig >= 0)
        for(int i = 0; i < numMigMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            migrationRate = migrationRate * (1 + sDriverMig * (1 - migrationRate / (maxRelMig + inputMigRate)) * rnd);
        }
    else
        for(int i = 0; i < numMigMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            migrationRate = migrationRate * (1 + sDriverMig * rnd);
        }
}
