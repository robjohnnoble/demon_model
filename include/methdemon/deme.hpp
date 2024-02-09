#ifndef DEME_HPP
#define DEME_HPP

#include "cell.hpp"
#include "parameters.hpp"

#include <string>
#include <vector>

class Deme {
private:
    // Fixed properties
    int K; // carrying capacity of the deme
    std::string side; // Left or right
    int identity; // Identity of the deme (index in tumour)
    float originTime;
    // Variable properties
    int population; // Number of cancer cells in the deme
    std::vector<Cell> cellList; // List of cells in the deme
    std::vector<float> avgMethArray; // Average methylation array of the deme
    int fissions; // fissions since the initial deme
    // rates
    float deathRate; // Death rate of cells in the deme (population dependent)
    float sumBirthRates; // Sum of birth rates of the cells in the deme
    float sumMigRates; // Sum of migration rates of the cells in the deme
    const float baseDeathRate; // Base death rate of cells in the deme (from input params)
public:
    // Constructor
    Deme(int K, std::string side, int identity, int population, int fissions, float deathRate, float baseDeathRate, float sumBirthRates, float sumMigrationRates);
    // Initialise first deme
    void initialise(std::shared_ptr<Genotype> firstGenotype, const InputParameters& params, const DerivedParameters& d_params);
    // Deme property handling
    void increment(int increment);
    void calculateAverageArray();
    // Deme events
    Deme demeFission(float originTime, bool firstFission=false);
    void pseudoFission();
    void moveCells(Deme& targetDeme);
    // Cell events
    int chooseCell();
    void cellDivision(int parentIndex, int* next_cell_id, int* nextGenotypeID, float gensElapsed, const InputParameters& params);
    void cellDeath(int cellIndex);
    // Rates handling
    void calculateSumsOfRates();
    // Getters
    int getK() const { return K; }
    std::string getSide() const { return side; }
    int getPopulation() const { return population; }
    int getIdentity() const { return identity; }
    float getDeathRate() const { return deathRate; }
    float getSumBirthRates() const { return sumBirthRates; }
    float getSumMigrationRates() const { return sumMigRates; }
    float getSumOfRates() const { return sumBirthRates + sumMigRates + population * deathRate; }
    float getCellBirth(int chosenCell) const { return cellList[chosenCell].getBirthRate(); }
    float getCellMig(int chosenCell) const { return cellList[chosenCell].getMigrationRate(); }
    float getOriginTime() const { return originTime; }
    int getFissions() const { return fissions; }
    std::vector<float> getAverageArray() const { return avgMethArray; }
    // Setters
    void setSide(std::string side) { this->side = side; }
    void setDeathRate() { this->deathRate = population > K ? baseDeathRate + 10 : baseDeathRate; }
    void setOriginTime(float originTime) { this->originTime = originTime; }
};

#endif // DEME_HPP
