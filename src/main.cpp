#include "Input.hpp"
#include "Initialise.hpp"
#include "RunSim.hpp"

int main(int argc, char *argv[]) {
    std::string input_and_output_path = get_input_path(argc, argv);
    std::cout << "Input and output directory: " << input_and_output_path << std::endl;
    std::string config_file_with_path = input_and_output_path + argv[2];
    std::cout << "Config file path: " << config_file_with_path << std::endl;

    boost::property_tree::ptree pt;
    boost::property_tree::info_parser::read_info(config_file_with_path, pt);

    InputParameters params = read_parameters(pt, config_file_with_path);
    run_sim(input_and_output_path, config_file_with_path, params);

    return 0;
}