#include "deme.hpp"

/////// Constructor
Deme::Deme(int K, std::string side, int identity, int population,
    int fissions, float deathRate, float sumBirthRates, float sumMigRates)
    : K(K), side(side), identity(identity), population(population), fissions(fissions),
        deathRate(deathRate), sumBirthRates(sumBirthRates), sumMigRates(sumMigRates) {
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
// increment or decrement population of the deme and update death rate
void Deme::increment(int increment, const InputParameters& params, std::string origin) {
    population += increment;
    setDeathRate(params);

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
// cell division
void Deme::cellDivision(int parentIndex, int* next_cell_id, const InputParameters& params, RandomNumberGenerator& rng) {
    Cell& parent = cellList[parentIndex];
    Cell daughter = Cell((*next_cell_id)++, parent.getGenotype(), parent.getDriverIndex(), identity, parent.getNumMeth(), parent.getNumDemeth(),
        parent.getFCpGs(), parent.getMethArray());
    parent.methylation(params);
    daughter.methylation(params);
    parent.mutation();
    daughter.mutation();
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
