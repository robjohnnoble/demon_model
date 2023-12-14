#include "cell.hpp"

/////// Constructor
Cell::Cell(int identity, std::shared_ptr<Genotype> genotype, int deme, int numMeth, int numDemeth, int fcpgs, std::vector<int> methArray)
    : identity(identity), genotype(genotype), deme(deme), numMeth(numMeth), numDemeth(numDemeth), fcpgs(fcpgs), methArray(methArray) {}

/////// Methylation array handling
// generate initial methylation array
void Cell::initialArray(const float manualArray) {
    methArray = std::vector<int>(fcpgs, 0);
    if(manualArray == -1) {
        for (int i = 0; i < fcpgs; i++) {
            double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
            rnd > 0.5 ? methArray[i] = 1 : methArray [i] = 0;
        }
    }
    else {
        for (int i = 0; i < fcpgs * manualArray; i++) {
            methArray[i] = 0;
        }
        for (int i = std::ceil(fcpgs * manualArray); i < fcpgs; i++) {
            double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
            rnd > 0.5 ? methArray[i] = 1 : methArray [i] = 0;
        }
    }
}
// methylation event
void Cell::methylation(const InputParameters& params) {
    for (int i = 0; i < fcpgs; i++) {
        double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
        int condition1 = methArray[i] == 0 && rnd < params.meth_rate;
        int condition2 = methArray[i] == 1 && rnd < params.demeth_rate;

        methArray[i] = methArray[i] + condition1 - condition2;
        numMeth += condition1;
        numDemeth += condition2;
    }
}

/////// Mutations
// mutation event
void Cell::mutation(int* next_genotype_id, const InputParameters& params) {
    int newBirthMut = RandomNumberGenerator::getInstance().poissonDist(genotype->getMuDriverBirth());
    int newMigMut = RandomNumberGenerator::getInstance().poissonDist(genotype->getMuDriverMig());

    if (newBirthMut || newMigMut) {
        std::shared_ptr<Genotype> newGenotype = std::make_shared<Genotype>(genotype->getIdentity(), (*next_genotype_id)++, genotype->getNumBirthMut() + newBirthMut, genotype->getNumMigMut() + newMigMut,
            genotype->getBirthRate(), genotype->getMigrationRate(), genotype->getMuDriverBirth(), genotype->getMuDriverMig(), genotype->getImmortal(), genotype->getOriginTime());
        newGenotype->setBirthRate(params);
        newGenotype->setMigrationRate(params);
        genotype = newGenotype;
    }
}