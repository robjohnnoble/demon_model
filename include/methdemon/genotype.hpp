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
    bool immortal; // indicator whether genotype record can be overwritten
    float muDriverBirth; // birth driver mutation rate
    float muDriverMig; // migration driver mutation rate
    // Time
    float originTime; // generation in which the genotype was created
public:
    // Constructor
    Genotype(int parent, int identity, int numBirthMut, int numMigMut, float birthRate, float migrationRate, float muDriverBirth, float muDriverMig, bool immortal, float originTime);
    // Cell birth events
    void increment(int increment);
    void addCell(int cellIdentity);
    std::vector<int> newMutations();
    void removeCell(int cellIdentity);
    // Rates handling
    void setBirthRate(const InputParameters& params);
    void setMigrationRate(const InputParameters& params);
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
    bool getImmortal() const { return immortal; }
    float getOriginTime() const { return originTime; }
};

#endif // GENOTYPE_HPP