#include "deme.hpp"

/////// Constructor
Deme::Deme(int K, std::string side, int identity, int population, int fissions, float deathRate, float sumBirthRates, float sumMigRates) : K(K), side(side), identity(identity), population(population), fissions(fissions), deathRate(deathRate), sumBirthRates(sumBirthRates), sumMigRates(sumMigRates) {
    avgMethArray.clear();
}

/////// Initialise first deme
void Deme::initialise(const InputParameters& params, const DerivedParameters& d_params) {
    // initialise first cell
    std::vector<int> tmpArray(d_params.fcpgs, 0);
    std::shared_ptr<Genotype> firstGenotype = std::make_shared<Genotype>(0, 0, 0, 0, 1, params.init_migration_rate, params.mu_driver_birth, params.mu_driver_migration, true, 0);
    Cell firstCell = Cell(0, firstGenotype, identity, 0, 0, d_params.fcpgs, tmpArray);
    firstCell.initialArray(params.manual_array);
}

/////// Deme property handling
// increment or decrement population of the deme and update all rates
void Deme::increment(int increment, const InputParameters& params, std::string origin) {
    population += increment;
    setDeathRate(params);
    calculateSumsOfRates();

    // check population sum
    int num_clones_in_deme = cellList.size();
    if (population != num_clones_in_deme) {
        std::cout << "ERROR: Population does not equal number of clones in deme." << std::endl
        << "deme identity: " << identity 
        << "; population: " << population << "; num_clones_in_deme: " << num_clones_in_deme << std::endl
        << "; increment: " << increment 
        << "; call origin: " << origin << std::endl;
        exit(1);
    }
}
// calculate the average methylation array of the deme
void Deme::calculateAverageArray(int fcpgs) {
    avgMethArray = std::vector<float>(fcpgs, 0);
    for (int i = 0; i < population; i++) {
        for (int j = 0; j < fcpgs; j++) {
            avgMethArray[j] += cellList[i].getFCpGSite(j);
        }
    }
    for (int i = 0; i < fcpgs; i++) {
        avgMethArray[i] /= population;
    }
}

/////// Cell events
// choose cell for event
int Deme::chooseCell() {
    std::vector<double> cumRates;
    double r;

    if (population == 1) {
        return 0;
    } else {
        double rateSum = 0.0;
        for (int i = 0; i < population; i++) {
            rateSum += deathRate + cellList[i].getBirthRate() + cellList[i].getMigrationRate();
            cumRates.push_back(rateSum);
        }
        double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
        r = rnd * cumRates.back();
    }

    if (population == 2) {
        return r < cumRates[0] ? 0 : 1;
    } else {
        auto it = std::lower_bound(cumRates.begin(), cumRates.end(), r);
        return std::distance(cumRates.begin(), it);
    }
}
// cell division
void Deme::cellDivision(int parentIndex, int* nextCellID, int* nextGenotypeID, const InputParameters& params, RandomNumberGenerator& rng) {
    Cell& parent = cellList[parentIndex];
    Cell daughter = Cell((*nextCellID)++, parent.getGenotype(), identity, parent.getNumMeth(), parent.getNumDemeth(), parent.getFCpGs(), parent.getMethArray());
    parent.methylation(params);
    daughter.methylation(params);
    parent.mutation(nextGenotypeID, params);
    daughter.mutation(nextGenotypeID, params);
    cellList.push_back(daughter);
}
// cell death
void Deme::cellDeath(int cellIndex) {
    std::swap(cellList[cellIndex], cellList.back());
    cellList.pop_back();
}

/////// Rates handling
// set death rate per cell based on current population
void Deme::setDeathRate(const InputParameters& params) {
    if (population <= K) {
        deathRate = params.baseline_death_rate;
    }
    else {
        deathRate = params.baseline_death_rate + 100;
    }
}
// calculate all rates
void Deme::calculateSumsOfRates() {
    sumBirthRates = 0;
    sumMigRates = 0;
    for (int i = 0; i < population; i++) {
        sumBirthRates += cellList[i].getBirthRate();
        sumMigRates += cellList[i].getMigrationRate();
    }
}