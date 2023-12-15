#ifndef CELL_HPP
#define CELL_HPP

#include "distributions.hpp"
#include "genotype.hpp"
#include "parameters.hpp"
#include <vector>

class Cell {
private:
    // Properties
    int identity; // Identity of the cell
    std::shared_ptr<Genotype> genotype; // Driver genotype of the cell
    // int driverIndex; // Index of the driver genotype in the tumour
    int deme; // Index of the deme in which the cell is located
    // numbers of methylation and demethylation events since the initial array
    int numMeth; // number of methylation events since initial array
    int numDemeth; // number of demethylation events since initial array
    // methylation array
    int fcpgs; // number of fCpG sites per cell
    std::vector<int> methArray; // fCpG array of the cell
    const float methRate; // methylation rate
    const float demethRate; // demethylation rate
public:
    // Constructor
    Cell(int identity, std::shared_ptr<Genotype> genotype, int deme, int numMeth, int numDemeth, int fcpgs, std::vector<int> methArray, float methRate, float demethRate);
    // Move constructor
    Cell(Cell&& other) noexcept : identity(other.identity), genotype(other.genotype), deme(other.deme), numMeth(other.numMeth), numDemeth(other.numDemeth), fcpgs(other.fcpgs), methArray(other.methArray), methRate(other.methRate), demethRate(other.demethRate) {}
    // Move assignment operator
    Cell& operator=(Cell&& other) noexcept;
    // Copy constructor
    Cell(const Cell& other);
    // Copy assignment operator
    Cell& operator=(const Cell& other);
    // Methylation array handling
    void initialArray(const float manualArray);
    void methylation();
    // Mutations
    void mutation(int* next_genotype_id, float gensElapsed, const InputParameters& params);
    // Getters
    int getIdentity() const { return identity; }
    std::shared_ptr<Genotype> getGenotype() const { return genotype; }
    // int getDriverIndex() const { return driverIndex; }
    int getDeme() const { return deme; }
    int getNumMeth() const { return numMeth; }
    int getNumDemeth() const { return numDemeth; }
    int getFCpGSite(int j) const { return methArray[j]; }
    int getFCpGs() const { return fcpgs; }
    float getMethRate() const { return methRate; }
    float getDemethRate() const { return demethRate; }
    float getBirthRate() const { return genotype->getBirthRate(); }
    float getMigrationRate() const { return genotype->getMigrationRate(); }
    std::vector<int> getMethArray() const { return methArray; }
    // Setters
    void setDeme(int deme) { this->deme = deme; }
};

#endif // CELL_HPP