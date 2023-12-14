#ifndef DEME_HPP
#define DEME_HPP

#include "parameters.hpp"
#include "distributions.hpp"
#include "cell.hpp"
#include <iostream>
#include <string>
#include <vector>

class Deme {
private:
    // Fixed properties
    int K; // carrying capacity of the deme
    std::string side; // Left or right
    int identity; // Identity of the deme (index in tumour)
    // Variable properties
    int population; // Number of cancer cells in the deme
    std::vector<Cell> cellList; // List of cells in the deme    
    std::vector<float> avgMethArray; // Average methylation array of the deme
    int fissions; // fissions since the initial deme
    // rates
    float deathRate; // Death rate of cells in the deme (population dependent)
    float sumBirthRates; // Sum of birth rates of the cells in the deme
    float sumMigRates; // Sum of migration rates of the cells in the deme
public:
    // Constructor
    Deme(int K, std::string side, int identity, int population, int fissions, float deathRate, float sumBirthRates, float sumMigrationRates);
    // Initialise first deme
    void initialise(const InputParameters& params, const DerivedParameters& d_params);
    // Deme property handling
    void increment(int increment, const InputParameters& params, std::string origin);
    void calculateAverageArray(int fcpgs);
    // Cell events
    void cellDivision(int parentIndex, int* next_cell_id, const InputParameters& params, RandomNumberGenerator& rng);
    void cellDeath(int cellIndex);
    // Rates handling
    void setDeathRate(const InputParameters& params);
    // Getters
    int getPopulation() const { return population; }
    int getIdentity() const { return identity; }
    float getDeathRate() const { return deathRate; }
    float getSumBirthRates() const { return sumBirthRates; }
    float getSumMigrationRates() const { return sumMigRates; }
    float getSumOfRates() const { return sumBirthRates + sumMigRates + population * deathRate; }
};

#endif // DEME_HPP