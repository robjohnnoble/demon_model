#include "Input.hpp"

// get input and path from terminal
std::string get_input_path(int argc, char *argv[]) {
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
InputParameters read_parameters(const boost::property_tree::ptree& pt, const std::string& config_file_path) {
    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(config_file_path, pt);

    InputParameters params;

    params.log2_deme_carrying_capacity = pt.get<int>("capacity.log2_deme_carrying_capacity");

    params.init_migration_rate = pt.get<float>("dispersal.init_migration_rate");
    params.migration_rate_scales_with_K = pt.get<int>("dispersal.migration_rate_scales_with_K");

    params.time0 = pt.get<int>("fission_times.time0");
    params.time1 = pt.get<int>("fission_times.time1");
    params.time2 = pt.get<int>("fission_times.time2");
    params.time3 = pt.get<int>("fission_times.time3");
    params.time4 = pt.get<int>("fission_times.time4");
    params.time5 = pt.get<int>("fission_times.time5");
    params.time6 = pt.get<int>("fission_times.time6");

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

    params.init_pop = pt.get<int>("initial_conditions.init_pop");

    params.record_matrix = pt.get<int>("output_indicators.record_matrix");
    params.write_demes_file = pt.get<int>("output_indicators.write_demes_file");
    params.write_clones_file = pt.get<int>("output_indicators.write_clones_file");
    params.write_phylo = pt.get<int>("output_indicators.write_phylo");
    params.calculate_total_diversity = pt.get<int>("output_indicators.calculate_total_diversity");

    return params;
}