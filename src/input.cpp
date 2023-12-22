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

    params.log2_deme_carrying_capacity = pt.get<int>("capacity.log2_deme_carrying_capacity");

    params.init_migration_rate = pt.get<float>("dispersal.init_migration_rate");
    params.migration_rate_scales_with_K = pt.get<int>("dispersal.migration_rate_scales_with_K");

    params.t0 = pt.get<int>("fission_times.t0");
    params.tL1 = pt.get<int>("fission_times.tL1");
    params.tL2 = pt.get<int>("fission_times.tL2");
    params.tL3 = pt.get<int>("fission_times.tL3");
    params.tR1 = pt.get<int>("fission_times.tR1");
    params.tR2 = pt.get<int>("fission_times.tR2");
    params.tR3 = pt.get<int>("fission_times.tR3");

    params.normal_birth_rate = pt.get<float>("fitness.normal_birth_rate");
    params.baseline_death_rate = pt.get<float>("fitness.baseline_death_rate");
    params.s_driver_birth = pt.get<float>("fitness.s_driver_birth");
    params.s_driver_migration = pt.get<float>("fitness.s_driver_migration");
    params.max_relative_birth_rate = pt.get<float>("fitness.max_relative_birth_rate");
    params.max_relative_migration_rate = pt.get<float>("fitness.max_relative_migration_rate");

    params.mu_driver_birth = pt.get<float>("mutation.mu_driver_birth");
    params.mu_driver_migration = pt.get<float>("mutation.mu_driver_migration");

    params.meth_rate = pt.get<float>("methylation.meth_rate");
    params.demeth_rate = pt.get<float>("methylation.demeth_rate");
    params.fCpG_loci_per_cell = pt.get<int>("methylation.fCpG_loci_per_cell");
    params.manual_array = pt.get<float>("methylation.manual_array");

    params.seed = pt.get<int>("seed.seed");

    params.max_time = pt.get<int>("stopping_conditions.max_time");
    params.max_generations = pt.get<int>("stopping_conditions.max_generations");
    params.max_fissions = pt.get<int>("stopping_conditions.max_fissions");

    params.init_pop = pt.get<int>("initial_conditions.init_pop");

    params.write_demes_file = pt.get<int>("output_indicators.write_demes_file");
    params.write_clones_file = pt.get<int>("output_indicators.write_clones_file");

    return params;
}