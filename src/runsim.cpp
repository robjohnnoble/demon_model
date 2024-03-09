#include "runsim.hpp"

float calculateTime(Tumour& tumour) {
    // implement calculations for gensAdded
    float tmp = tumour.sumAllRates();
    float rnd = RandomNumberGenerator::getInstance().expDist(1);

    return rnd / tmp;
}

void runSim(const std::string& input_and_output_path,
    const std::string& config_file_with_path, const InputParameters& params) {
    int iterations = 0;
    float outputTimer = 0;
    float gensAdded; // time tracking
    // derive derived parameters
    DerivedParameters d_params = deriveParameters(params);
    // initialise output files
    FileOutput finalDemes(input_and_output_path + "final_demes.csv");
    finalDemes.writeDemesHeader();
    // initialise tumour
    Tumour tumour(params, d_params);
    // NOTE: Implement event counter eventually (not that important tbh)
    std::cout << "Initialised simulation." << std::endl;
    // start timer
    auto start = std::chrono::high_resolution_clock::now();
    // sort out column headers in output files
    while(tumour.getFissionsPerDeme() < params.max_fissions ||
            tumour.getNumDemes() < d_params.max_demes) {
        tumour.event(params, d_params);

        // update time
        iterations++;
        gensAdded = calculateTime(tumour);
        tumour.setGensElapsed(gensAdded);
        outputTimer += gensAdded;

        // write to stdout and files every 5 generations
        if(outputTimer >= 5) {
            int numGenotypes = tumour.getNumGenotypes();
            int numDemes = tumour.getNumDemes();
            int numCells = tumour.getNumCells();
            std::cout << "Generations elapsed: " << tumour.getGensElapsed()
                      << ", Iterations: " << iterations
                      << ", Mean Fissions: " << tumour.getFissionsPerDeme()
                      << std::endl;
            std::cout << "Number of cells: " << numCells << std::endl;
            std::cout << "Number of driver genotypes: " << numGenotypes
                      << std::endl;
            std::cout << "Number of demes: " << numDemes << std::endl;
            std::cout << tumour.getNextCellID() << " cells ever created; "
                      << tumour.getNextGenotypeID()
                      << " genotypes ever created." << std::endl;
            outputTimer = 0;
            finalDemes.writeDemesFile(tumour);
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "End of simulation." << std::endl
    << tumour.getNumDemes() << " demes; " << tumour.getNumCells() << " cells; "
    << tumour.getGensElapsed() << " generations; "
    << tumour.getFissionsPerDeme() << " mean fissions per deme." << std::endl;
    std::cout << "Running time: " << elapsed.count() << " seconds." << std::endl;
    finalDemes.writeDemesFile(tumour);
}
