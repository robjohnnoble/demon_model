#include "cell.hpp"

/////// Constructor
Cell::Cell(int identity, std::shared_ptr<Genotype> genotype, int deme, int numMeth, int numDemeth, int fcpgs, std::vector<int> methArray, float methRate, float demethRate)
    : identity(identity), genotype(genotype), deme(deme), numMeth(numMeth), numDemeth(numDemeth), fcpgs(fcpgs), methArray(methArray), methRate(methRate), demethRate(demethRate) {}
// Move assignment operator
Cell& Cell::operator=(Cell&& other) noexcept {
    // Guard against self-assignment
    if (this != &other) {
        // Transfer ownership of 'other's resources to 'this'
        identity = other.identity;
        genotype = std::move(other.genotype);
        deme = other.deme;
        numMeth = other.numMeth;
        numDemeth = other.numDemeth;
        fcpgs = other.fcpgs;
        methArray = std::move(other.methArray);
    }
    return *this;
}
// Copy constructor
Cell::Cell(const Cell& other)
    : identity(other.identity),
      genotype(other.genotype),
      deme(other.deme),
      numMeth(other.numMeth),
      numDemeth(other.numDemeth),
      fcpgs(other.fcpgs),
      methArray(other.methArray),
      methRate(other.methRate),
      demethRate(other.demethRate) {}
// Copy assignment operator
Cell& Cell::operator=(const Cell& other) {
    if (this != &other) { // Guard against self-assignment
        identity = other.identity;
        genotype = other.genotype; // Assuming shared_ptr should be copied
        deme = other.deme;
        numMeth = other.numMeth;
        numDemeth = other.numDemeth;
        fcpgs = other.fcpgs;
        methArray = other.methArray; // std::vector supports direct copy
        // Note: No need to assign methRate and demethRate as they are const
    }
    return *this;
}

/////// Methylation array handling
// generate initial methylation array
void Cell::initialArray(const float manualArray) {
    methArray = std::vector<int>(fcpgs, 0);
    if(manualArray == -1) {
        for (int i = 0; i < fcpgs; i++) {
            double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
            rnd > 0.5 ? methArray[i] = 1 : methArray[i] = 0;
        }
    }
    else {
        for (int i = 0; i < fcpgs * manualArray; i++) {
            methArray[i] = 0;
        }
        for (int i = std::ceil(fcpgs * manualArray); i < fcpgs; i++) {
            double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
            rnd > 0.5 ? methArray[i] = 1 : methArray[i] = 0;
        }
    }
    // for (int i = 0; i < fcpgs; i++) {
    //     std::cout << methArray[i];
    // }
    // std::cout << std::endl;
}
// methylation event
void Cell::methylation() {
    for (int i = 0; i < fcpgs; i++) {
        double rnd = RandomNumberGenerator::getInstance().unitUnifDist();
        int condition1 = methArray[i] == 0 && rnd < methRate;
        int condition2 = methArray[i] == 1 && rnd < demethRate;

        methArray[i] = methArray[i] + condition1 - condition2;
        numMeth += condition1;
        numDemeth += condition2;
    }
    // std::cout << "methylations: " << numMeth << "; demethylations: " << numDemeth << std::endl;
}

/////// Mutations
// mutation event
void Cell::mutation(int *next_genotype_id, float gensElapsed,
                    const InputParameters &params) {
  int newBirthMut = RandomNumberGenerator::getInstance().poissonDist(
      genotype->getMuDriverBirth());
  int newMigMut = RandomNumberGenerator::getInstance().poissonDist(
      genotype->getMuDriverMig());

  if (newBirthMut || newMigMut) {
    std::shared_ptr<Genotype> newGenotype = std::make_shared<Genotype>(
        genotype->getIdentity(), (*next_genotype_id)++,
        genotype->getNumBirthMut() + newBirthMut,
        genotype->getNumMigMut() + newMigMut, 0, 0, gensElapsed, params);
    newGenotype->setBirthRate();
    newGenotype->setMigrationRate();
    genotype = newGenotype;
    }
}
