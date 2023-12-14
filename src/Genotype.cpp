#include "genotype.hpp"

/////// Constructor
Genotype::Genotype(int parent, int identity,
    int numBirthMut, int numMigMut, float birthRate, float migrationRate,
    float muDriverBirth, float muDriverMig, bool immortal, float originTime)
    : parent(parent), identity(identity),
    numBirthMut(numBirthMut), numMigMut(numMigMut), birthRate(birthRate), migrationRate(migrationRate),
    muDriverBirth(muDriverBirth), muDriverMig(muDriverMig), immortal(immortal), originTime(originTime) {}

/////// Cell birh events
// increment population
// void Genotype::increment(int increment) {
//     population += increment;
//     if (population == 0) {
//         immortal = false;
//     }
// }
// // add cell to cell list
// void Genotype::addCell(int cell) {
//     cellList.push_back(cell);
// }
// // mutations
// std::vector<int> Genotype::newMutations() {
//     // draw the number of birth and migration drivers from Poisson distributions
//     int newBirth = RandomNumberGenerator::getInstance().poissonDist(muDriverBirth);
//     int newMig = RandomNumberGenerator::getInstance().poissonDist(muDriverMig);

//     // return the number of birth and migration drivers
//     std::vector<int> newDrivers = {newBirth, newMig};
//     return newDrivers;
// }

/////// Rates handling
// birth rate
void Genotype::setBirthRate(const InputParameters& params) {
    float advantage = params.s_driver_birth;
    birthRate = 1;

    if (params.max_relative_birth_rate >= 0) 
        for(int i = 0; i < numBirthMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            birthRate = birthRate * (1 + advantage * (1 - birthRate / params.max_relative_birth_rate) * rnd);
        }
    else
        for(int i = 0; i < numBirthMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            birthRate = birthRate * (1 + advantage * rnd);
        }
    
    if (birthRate >= params.max_relative_birth_rate + 100) {
        std::cout << "ERROR: Birth rate at carrying capacity exceeds death rate." << std::endl;
        exit(1);
    }
}

void Genotype::setMigrationRate(const InputParameters& params) {
    float advantage = params.s_driver_migration;
    migrationRate = params.init_migration_rate;

    if (params.max_relative_migration_rate >= 0) 
        for(int i = 0; i < numMigMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            migrationRate = migrationRate * (1 + advantage * (1 - migrationRate / (params.max_relative_migration_rate + params.init_migration_rate)) * rnd);
        }
    else
        for(int i = 0; i < numMigMut; i++) {
            float rnd = RandomNumberGenerator::getInstance().expDist(1);
            migrationRate = migrationRate * (1 + advantage * rnd);
        }
}