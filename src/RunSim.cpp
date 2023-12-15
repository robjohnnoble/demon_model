#include "runsim.hpp"
#include "output.hpp"

float calculateTime(Tumour& tumour) {
    // implement calculations for gensAdded
    float tmp = tumour.sumAllRates();
    float rnd = RandomNumberGenerator::getInstance().expDist(1);

    return rnd / tmp;
}

void runSim(const std::string& input_and_output_path,
    const std::string& config_file_with_path, const InputParameters& params) {
    int iterations;
    int timeTracker = 0;
    float gensElapsed = 0, gensAdded; // time tracking
    // derive derived parameters
    DerivedParameters d_params = deriveParameters(params);
    // initialise output files
    FileOutput finalDemes(input_and_output_path + "final_demes.txt");
    // initialise tumour
    Tumour tumour(params, d_params);
    // initialise event counter
    EventCounter event_counter;
    std::cout << "Initialised simulation." << std::endl;
    // open output files
    // sort out column headers in output files
    while(gensElapsed < params.max_generations &&
        tumour.fissionsPerDeme() < params.max_fissions) {
        // update all rate sums
        // tumour.calculate_all_rates(params, d_params);
        // choose deme, cell, event
        tumour.event(params);
        // perform event
        if(iterations % 10000 == 0) {
            int numGenotypes = tumour.getNumGenotypes();
            int numDemes = tumour.getNumDemes();
            int numCells = tumour.getNumCells();
            std::cout << "Generations elapsed: " << gensElapsed << std::endl;
            std::cout << "Number of cells: " << numCells << std::endl;
            std::cout << "Number of driver genotypes: " << numGenotypes << std::endl;
            std::cout << "Number of demes: " << numDemes << std::endl;
            std::cout << tumour.getNextCellID() << " cells ever created; " << tumour.getNextGenotypeID() << " genotypes ever created." << std::endl;
            timeTracker = 0;
        }
        // update time
        iterations++;
        gensAdded = calculateTime(tumour);
        gensElapsed += gensAdded;
        timeTracker += gensAdded;
    }


    std::cout << "End of simulation." << std::endl
    << tumour.getNumDemes() << " demes; " << tumour.getNumCells() << " cells; " 
    << tumour.fissionsPerDeme() << " mean fissions per deme." << std::endl;
    finalDemes.write(tumour);
}