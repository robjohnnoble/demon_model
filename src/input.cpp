#include "input.hpp"

// get input and path from terminal
std::string getInputPath(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Missing argument providing path to the config file and/or the config file name." << std::endl;
        exit(0);
    }

    std::string path = argv[1];
    if (path.back() != '/') {
        path += '/';
    }

    return path;
}

// read parameters from config file
InputParameters readParameters(const boost::property_tree::ptree& pt, const std::string& config_file_path) {
    InputParameters params;

    params.deme_carrying_capacity = pt.get<int>("capacity.deme_carrying_capacity");

    params.init_migration_rate = pt.get<float>("dispersal.init_migration_rate");
    params.left_demes = pt.get<int>("dispersal.left_demes");
    params.right_demes = pt.get<int>("dispersal.right_demes");
    params.migration_rate_scales_with_K = pt.get<int>("dispersal.migration_rate_scales_with_K");

    params.mu_driver_birth = pt.get<float>("mutation.mu_driver_birth");
    params.mu_driver_migration = pt.get<float>("mutation.mu_driver_migration");

    params.normal_birth_rate = pt.get<float>("fitness.normal_birth_rate");
    params.baseline_death_rate = pt.get<float>("fitness.baseline_death_rate");
    params.s_driver_birth = pt.get<float>("fitness.s_driver_birth");
    params.s_driver_migration = pt.get<float>("fitness.s_driver_migration");
    params.max_relative_birth_rate = pt.get<float>("fitness.max_relative_birth_rate");
    params.max_relative_migration_rate = pt.get<float>("fitness.max_relative_migration_rate");

    params.meth_rate = pt.get<float>("methylation.meth_rate");
    params.demeth_rate = pt.get<float>("methylation.demeth_rate");
    params.fCpG_loci_per_cell = pt.get<int>("methylation.fCpG_loci_per_cell");
    params.manual_array = pt.get<float>("methylation.manual_array");

    params.seed = pt.get<int>("rng_seed.seed");

    params.max_time = pt.get<int>("stopping_conditions.max_time");
    params.max_generations = pt.get<int>("stopping_conditions.max_generations");
    params.max_fissions = pt.get<int>("stopping_conditions.max_fissions");

    params.init_pop = pt.get<int>("initial_conditions.init_pop");
    params.fission_config = pt.get<int>("initial_conditions.fission_config");

    params.write_demes_file = pt.get<int>("output_indicators.write_demes_file");
    params.write_clones_file = pt.get<int>("output_indicators.write_clones_file");

    return params;
}
