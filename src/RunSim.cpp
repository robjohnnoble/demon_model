#include "RunSim.hpp"

void run_sim(const std::string& input_and_output_path,
    const std::string& config_file_with_path, const InputParameters& params) {
    // initialise random number generator
    RandomNumberGenerator rng;
    rng.set_seed(params.seed);

    // initialise derived parameters
    DerivedParameters d_params = derived_parameters(params);

    // initialise tumour
    Tumour tumour;
    tumour.initialise(params, d_params, rng);

    // initialise event counter
    EventCounter event_counter;
    
    // open output files
    // sort out column headers in output files

    while(tumour.check_time() < params.max_generations) {
        tumour.update_time(rng, params);

        // choose deme, weights determined by number of rates
        int chosen_deme = tumour.choose_deme(rng);

        // select clone within deme
        int chosen_clone = tumour.choose_clone(chosen_deme, rng);
        // select event type
        std::string event_type = tumour.choose_event_type(chosen_deme, chosen_clone, rng);

        // perform event
        if (event_type == "birth") {
            tumour.cell_division(event_counter, rng, chosen_deme, chosen_clone, params);
        }
        else if (event_type == "death") {
            tumour.cell_death(event_counter, chosen_deme, chosen_clone, params);
        }
        else if (event_type == "fission") {
            tumour.deme_fission(chosen_deme, event_counter, rng, params);
        }

        std::cout << "Event_type: " << event_type << "\n";
        std::cout << "Time: " << tumour.check_time() << "\n";
        std::cout << "Number of cells: " << tumour.num_cells() << "\n";
        std::cout << "Number of clones: " << tumour.num_clones() << "\n";
        std::cout << "Number of driver genotypes: " << tumour.num_driver_genotypes() << "\n";
        std::cout << "Number of demes: " << tumour.num_demes() << "\n";

        tumour.iterations++;
    }

    tumour.final_output();
}