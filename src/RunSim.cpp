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

        // choose deme, weights determined by number of rates
        int chosen_deme = tumour.choose_deme(rng);

        // select clone within deme
        int chosen_clone = tumour.choose_clone(chosen_deme, rng);
        // select event type
        std::string event_type = tumour.choose_event_type(chosen_deme, chosen_clone, rng);

        // perform event
        if (event_type == "birth") {
            tumour.cell_division(event_counter, rng, chosen_deme, chosen_clone, params);
            if(tumour.fission_ready(chosen_deme, rng, true)) {
                tumour.deme_fission(chosen_deme, event_counter, rng, params);
            }
        }
        else if (event_type == "death") {
            tumour.cell_death(event_counter, chosen_deme, chosen_clone, params);
        }
        else if (event_type == "fission" && tumour.fission_ready(chosen_deme, rng, false)) {
            tumour.deme_fission(chosen_deme, event_counter, rng, params);
        }

        if(tumour.iterations % 10000 == 0) {
            float time = tumour.check_time();
            int num_cells = tumour.num_cells();
            int num_clones = tumour.num_clones();
            int num_driver_genotypes = tumour.num_driver_genotypes();
            int num_demes = tumour.num_demes();
            std::cout << "Event_type: " << event_type << std::endl;
            std::cout << "Generations elapsed: " << time << std::endl;
            std::cout << "Number of cells: " << num_cells << std::endl;
            std::cout << "Number of clones: " << num_clones << std::endl;
            std::cout << "Number of driver genotypes: " << num_driver_genotypes << std::endl;
            std::cout << "Number of demes: " << num_demes << std::endl;
        }

        tumour.iterations++;
        tumour.update_time(rng, params);
    }

    //tumour.final_output();
    std::cout << "End of simulation." << std::endl;
}