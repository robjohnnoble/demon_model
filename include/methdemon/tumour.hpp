#ifndef TUMOUR_HPP
#define TUMOUR_HPP

#include "parameters.hpp"
#include "deme.hpp"
#include "genotype.hpp"
#include <vector>

class Tumour {
private:
    // cell containers
    std::vector<Deme> demes;
    std::vector<std::shared_ptr<Genotype> > genotypes;
    // cell and genotype ID tracking
    int nextGenotypeID = 1;
    int nextCellID = 1;
    // temporal variables
    float gensElapsed = 0;
    float outputTimer = 0;
    std::vector<float> fissionTimes;
    float* nextFissionL = 0;
    float* nextFissionR = 0;
    int maxGens = 0;
    // sums of rates
    // double sumBirthRates = 0;
    // double sumDeathRates = 0;
    // double sumMigRates = 0;
    // misc
    int fissionConfig = 0;
public:
    // Constructor
    Tumour(const InputParameters& params, const DerivedParameters& d_params);
    // // deme rates
    // void calculate_deme_birth_rate(Deme& deme);
    // void calculate_deme_migration_rate(Deme& deme);
    // void calculate_all_rates(const InputParameters& params, const DerivedParameters& d_params);
    // choose deme, cell and event type
    int chooseDeme();
    int chooseCell(int chosenDeme) { return demes[chosenDeme].chooseCell(); };
    std::string chooseEventType(int chosenDeme, int chosenCell);
    //perform event
    void event(const InputParameters& params);
    // sum all rates (for time tracking)
    float sumAllRates();
    // Getters
    int getNextCellID() const { return nextCellID; }
    int getNextGenotypeID() const { return nextGenotypeID; }
    int getNumCells() const;
    float fissionsPerDeme();
    int getNumDemes() const { return demes.size(); }
    int getNumGenotypes() const { return genotypes.size(); }
    float getGensElapsed() const { return gensElapsed; }
    float getOutputTimer() const { return outputTimer; }
    Deme& getDeme(int index) { return demes[index]; }
    // Setters
    void setGensElapsed(float gensAdded = 0) { gensElapsed += gensAdded; }
};

#endif // TUMOUR_HPP
