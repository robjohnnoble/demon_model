#include "tumour.hpp"

/////// Constructor
Tumour::Tumour(const InputParameters& params,
    const DerivedParameters& d_params) {
    demes.clear();
    genotypes.clear();
    // driver genotypes:
    std::shared_ptr<Genotype> firstGenotype= std::make_shared<Genotype>(0, 0, 0, 0, 1, params.init_migration_rate, 0, params);
    genotypes.push_back(firstGenotype);

    // demes:
    Deme firstDeme(d_params.K, "left", 0, 1, 0, params.baseline_death_rate, params.baseline_death_rate, 1, params.init_migration_rate);
    demes.push_back(firstDeme);
    demes.back().initialise(firstGenotype, params, d_params);

    // fission times
    fissionTimes.push_back(params.time0);
    fissionTimes.push_back(params.time1);
    fissionTimes.push_back(params.time2);
    fissionTimes.push_back(params.time3);
    fissionTimes.push_back(params.time4);
    fissionTimes.push_back(params.time5);
    fissionTimes.push_back(params.time6);

    nextFission = fissionTimes[0];
}

/////// Choose events based on rate sums
// choose deme
int Tumour::chooseDeme() {
    std::vector<double> cumRates(demes.size());
    double r;
    int res = 0;

    if (cumRates.size() == 1) {
        return res;
    }
    else {
        double sumOfRates = demes[0].getSumOfRates();
        cumRates[0] = sumOfRates;
        for (int i = 1; i < cumRates.size(); i++) {
            sumOfRates = demes[i].getSumOfRates();
            cumRates[i] += sumOfRates + cumRates[i - 1];
        }
        double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
        r = rnd * cumRates.back();
    }
    if (cumRates.size() == 2) {
        res = r < cumRates[0] ? 0 : 1;
        return res;
    }
    else {
        res = std::lower_bound(cumRates.begin(), cumRates.end(), r) - cumRates.begin();
        return res;
    }
}
// choose event type
std::string Tumour::chooseEventType(int chosenDeme, int chosenCell) {
    std::vector<float> cumRates;
    int ctr = 0;
    float res;
    double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
    // cell birth rate
    cumRates.push_back(demes[chosenDeme].getCellBirth(chosenCell));
    ctr++;
    // cell death rate
    cumRates.push_back(cumRates[ctr - 1] + demes[chosenDeme].getDeathRate());
    ctr++;
    // cell migration rate
    cumRates.push_back(cumRates[ctr - 1] + demes[chosenDeme].getCellMig(chosenCell));
    // weighted choice
    res = rnd * cumRates.back();
    if(res < cumRates[0]) {
        return "birth";
    }
    else if(res < cumRates[1]) {
        return "death";
    }
    else {
        return "fission";
    }
}
// perform event
void Tumour::event(const InputParameters& params) {
    int chosenDeme = chooseDeme();
    int chosenCell = demes[chosenDeme].chooseCell();
    std::string eventType = chooseEventType(chosenDeme, chosenCell);

    if (eventType == "birth") {
        demes[chosenDeme].cellDivision(chosenCell, &nextCellID, &nextGenotypeID, gensElapsed, params);
        float rnd = RandomNumberGenerator::getInstance().unitUnifDist();
        if (demes[chosenDeme].getPopulation() >= demes[chosenDeme].getK() && rnd < demes[chosenDeme].getSumMigrationRates())
            if (gensElapsed >= nextFission) {
                if (nextFission == fissionTimes[0]) {
                    nextFission = fissionTimes[1];
                    demes[chosenDeme].demeFission(true);
                }
                else {
                    demes[chosenDeme].demeFission();
                }
            }
            else {
                demes[chosenDeme].pseudoFission();
            }
    }
    else if (eventType == "death") {
        demes[chosenDeme].cellDeath(chosenCell);
    }
    else if (eventType == "fission") {
        if (demes[chosenDeme].getPopulation() >= demes[chosenDeme].getK())
            if (gensElapsed >= nextFission) {
                if (nextFission == fissionTimes[0]) {
                    nextFission = fissionTimes[1];
                    demes[chosenDeme].demeFission(true);
                }
                else {
                    demes[chosenDeme].demeFission();
                }
            }
            else {
                demes[chosenDeme].pseudoFission();
            }
        }
    else {
        std::cout << "Error: invalid event type" << std::endl;
    }
}

/////// Sum all rates
float Tumour::sumAllRates() {
    float res = 0;
    for (int i = 0; i < demes.size(); i++) {
        res += demes[i].getSumOfRates();
    }
    return res;
}

/////// Getters
// get total number of cells tracked in tumour
int Tumour::getNumCells() const {
    int res = 0;
    for (int i = 0; i < demes.size(); i++) {
        res += demes[i].getPopulation();
    }
    return res;
}