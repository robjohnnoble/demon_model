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
public:
    // Constructor
    Cell(int identity, std::shared_ptr<Genotype> genotype, int deme, int numMeth, int numDemeth, int fcpgs, std::vector<int> methArray);
    // Methylation array handling
    void initialArray(const float manualArray);
    void methylation(const InputParameters& params);
    // Mutations
    void mutation(int* next_genotype_id, const InputParameters& params);
    // Getters
    int getIdentity() const { return identity; }
    std::shared_ptr<Genotype> getGenotype() const { return genotype; }
    // int getDriverIndex() const { return driverIndex; }
    int getDeme() const { return deme; }
    int getNumMeth() const { return numMeth; }
    int getNumDemeth() const { return numDemeth; }
    int getFCpGSite(int j) const { return methArray[j]; }
    int getFCpGs() const { return fcpgs; }
    std::vector<int> getMethArray() const { return methArray; }
};

#endif // CELL_HPP