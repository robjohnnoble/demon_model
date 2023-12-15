#ifndef TUMOUR_HPP
#define TUMOUR_HPP

#include "parameters.hpp"
#include "distributions.hpp"
#include "deme.hpp"
#include "genotype.hpp"
#include "output.hpp"
#include "macros.hpp"
#include <set>
#include <functional>
#include <vector>
#include <iostream>

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
    std::vector<int> fissionTimes;
    int nextFission = 0;    
    // sums of rates
    double sumBirthRates = 0;
    double sumDeathRates = 0;
    double sumMigRates = 0;
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
    int getNumDemes() const { return demes.size(); }
    int getNumGenotypes() const { return genotypes.size(); }
    float getGensElapsed() const { return gensElapsed; }
    float getOutputTimer() const { return outputTimer; }
};

#endif // TUMOUR_HPP