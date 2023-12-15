#ifndef GENOTYPE_HPP
#define GENOTYPE_HPP

#include "distributions.hpp"
#include "parameters.hpp"
#include <iostream>

class Genotype {
private:
    // Properties
    int parent; // parent's unique ID
    int identity; // unique ID of the driver genotype
    // int index;
    // Cells
    // int population; // number of cells with this genotype
    // std::vector<int> cellList; // identities of cells with this genotype
    // Numbers of mutations and methylation events
    int numBirthMut; // number of driver mutations
    int numMigMut; // number of migration mutations
    // Rates
    float birthRate; // birth rate of the genotype
    float migrationRate; // migration rate of the genotype
    // float baseMigRate; // base migration rate
    const float muDriverBirth; // birth driver mutation rate
    const float muDriverMig; // migration driver mutation rate
    const float inputMigRate; // input migration rate
    const float maxRelBirth; // maximum relative birth rate
    const float maxRelMig; // maximum relative migration rate
    const float sDriverBirth; // advantage of birth driver mutations
    const float sDriverMig; // advantage of migration driver mutations
    // Time
    float originTime; // generation in which the genotype was created
public:
    // Constructor
    Genotype(int parent, int identity, int numBirthMut, int numMigMut, float birthRate, float migrationRate, float originTime, const InputParameters& params);
    // Cell birth events
    std::vector<int> newMutations();
    void removeCell(int cellIdentity);
    // Getters
    // int getPopulation() const { return population; }
    int getParent() const { return parent; }
    int getIdentity() const { return identity; }
    // int getIndex() const { return index; }
    int getNumBirthMut() const { return numBirthMut; }
    int getNumMigMut() const { return numMigMut; }
    float getBirthRate() const { return birthRate; }
    float getMigrationRate() const { return migrationRate; }
    float getMuDriverBirth() const { return muDriverBirth; }
    float getMuDriverMig() const { return muDriverMig; }
    float getOriginTime() const { return originTime; }
    // Setters
    void setBirthRate();
    void setMigrationRate();
    // void setBaseMigRate();
};

#endif // GENOTYPE_HPP