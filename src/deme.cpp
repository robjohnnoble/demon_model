#include "deme.hpp"

/////// Constructor
Deme::Deme(int K, std::string side, int identity, int population, int fissions, float deathRate, float baseDeathRate, float sumBirthRates, float sumMigRates) : K(K), side(side), identity(identity), population(population), fissions(fissions), deathRate(deathRate), sumBirthRates(sumBirthRates), sumMigRates(sumMigRates), baseDeathRate(baseDeathRate) {
    avgMethArray.clear();
    cellList.clear();
}

/////// Initialise first deme
void Deme::initialise(std::shared_ptr<Genotype> firstGenotype, const InputParameters& params, const DerivedParameters& d_params) {
    // initialise first cell
    std::vector<int> tmpArray(d_params.fcpgs, 0);
    Cell firstCell = Cell(0, firstGenotype, identity, 0, 0, d_params.fcpgs, tmpArray, params.meth_rate, params.demeth_rate);
    firstCell.initialArray(params.manual_array);
    cellList.push_back(std::move(firstCell));
    calculateAverageArray();
}

/////// Deme property handling
// increment or decrement population of the deme and update all rates
void Deme::increment(int increment) {
    population += increment;
    setDeathRate();
    calculateSumsOfRates();

    // check population sum
    int num_clones_in_deme = cellList.size();
    if (population != num_clones_in_deme) {
        std::cout << "ERROR: Population does not equal number of clones in deme." << std::endl
        << "deme identity: " << identity
        << "; population: " << population << "; num_clones_in_deme: " << num_clones_in_deme << std::endl
        << "; increment: " << increment << std::endl;
        // << "; call origin: " << origin << std::endl;
        exit(1);
    }
}
// calculate the average methylation array of the deme
void Deme::calculateAverageArray() {
    int fcpgs = cellList[0].getFCpGs();
    avgMethArray = std::vector<float>(fcpgs / 2, 0);
    for (int i = 0; i < population; i++) {
        for (int j = 0; j < fcpgs / 2; j++) {
            avgMethArray[j] += static_cast<float>(cellList[i].getFCpGSite(j) + cellList[i].getFCpGSite(j + fcpgs / 2)) / 2.0;
        }
    }
    for (int i = 0; i < fcpgs / 2; i++) {
        avgMethArray[i] /= static_cast<float>(population);
    }
}

/////// Deme events
// deme fission - returns new deme
Deme Deme::demeFission(float originTime, bool firstFission) {
    fissions++;
    // initialise new deme
    Deme newDeme = Deme(K, side, identity + 1, 0, fissions, 0, baseDeathRate, 0, 0);
    if (firstFission) newDeme.setSide("right");
    moveCells(newDeme);
    // update origin deme
    setDeathRate();
    calculateSumsOfRates();
    calculateAverageArray();
    // update new deme
    newDeme.setDeathRate();
    newDeme.calculateSumsOfRates();
    newDeme.calculateAverageArray();
    newDeme.setOriginTime(originTime);
    return newDeme;
}
// pseudo fission - kill half the population randomly
void Deme::pseudoFission() {
    fissions++;
    int numCellsToKill = RandomNumberGenerator::getInstance().stochasticRound(population / 2.0);
    // get set of indices of cells to kill
    std::vector<int> indices;
    for (int i = 0; i < population; i++) {
        indices.push_back(i);
    }
    // shuffle indices
    std::shuffle(indices.begin(), indices.end(), RandomNumberGenerator::getInstance().getEngine());
    // remove the first `numCellsToKill` indices in descending order
    std::sort(indices.begin(), indices.begin() + numCellsToKill, std::greater<int>());
    for (int i = 0; i < numCellsToKill; i++) {
        cellDeath(indices[i]);
    }
    setDeathRate();
    calculateSumsOfRates();
}
// move cells to target deme
void Deme::moveCells(Deme& targetDeme) {
    int numCellsToMove = RandomNumberGenerator::getInstance().stochasticRound(population / 2.0);
    // get set of indices of cells to move
    std::vector<int> indices;
    for (int i = 0; i < population; i++) {
        indices.push_back(i);
    }
    // shuffle indices
    std::shuffle(indices.begin(), indices.end(), RandomNumberGenerator::getInstance().getEngine());
    // move the first `numCellsToMove` indices in descending order
    std::sort(indices.begin(), indices.begin() + numCellsToMove, std::greater<int>());
    for (int i = 0; i < numCellsToMove; i++) {
        int index = indices[i];
        targetDeme.cellList.push_back(std::move(cellList[index]));
        targetDeme.cellList.back().setDeme(targetDeme.getIdentity());
        if(index != population - 1) {
            std::swap(cellList[index], cellList.back());
        }
        cellList.pop_back();
    }
    increment(-numCellsToMove);
    targetDeme.increment(numCellsToMove);
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
void Deme::cellDivision(int parentIndex, int* nextCellID, int* nextGenotypeID, float gensElapsed, const InputParameters& params) {
    Cell& parent = cellList[parentIndex];
    Cell daughter = Cell((*nextCellID)++, parent.getGenotype(), identity, parent.getNumMeth(), parent.getNumDemeth(), parent.getFCpGs(), parent.getMethArray(), parent.getMethRate(), parent.getDemethRate());
    parent.methylation();
    daughter.methylation();
    parent.mutation(nextGenotypeID, gensElapsed, params);
    daughter.mutation(nextGenotypeID, gensElapsed, params);
    cellList.push_back(std::move(daughter));
    increment(1);
}
// cell death
void Deme::cellDeath(int cellIndex) {
    std::swap(cellList[cellIndex], cellList.back());
    cellList.pop_back();
    increment(-1);
}

/////// Rates handling
// calculate all rates
void Deme::calculateSumsOfRates() {
    sumBirthRates = 0;
    sumMigRates = 0;
    for (int i = 0; i < population; i++) {
        sumBirthRates += cellList[i].getBirthRate();
        sumMigRates += cellList[i].getMigrationRate();
    }
}
