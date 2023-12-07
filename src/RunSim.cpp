#include "RunSim.hpp"

void run_sim(const std::string& input_and_output_path,
    const std::string config_file_with_path, const InputParameters& params) {
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
        Deme chosen_deme = tumour.choose_deme(rng);

        // select clone within deme
        Clone chosen_clone = tumour.choose_clone(chosen_deme, rng);
        // select event type
        event_type = tumour.choose_event_for_clone(chosen_deme, chosen_clone);

        // perform event
        if (event_type == "birth") {
            tumour.cell_division(chosen_clone, event_counter, rng, chosen_deme);
        }
        else if (event_type == "death") {
            tumour.cell_death(chosen_clone, event_counter, rng, chosen_deme);
        }
        else if (event_type == "fission") {
            tumour.deme_fission(chosen_clone, event_counter, rng, chosen_deme);
        }
        else {
            std::cout << "Error: event type not recognised" << std::endl;
        }
    }

    tumour.final_output();
}