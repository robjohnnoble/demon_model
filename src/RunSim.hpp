#ifndef RUNSIM_HPP
#define RUNSIM_HPP

#include "Initialise.hpp"
#include "Tumour.hpp"

#include <string>
#include <vector>

void run_sim(const std::string& input_and_output_path, const std::string& config_file_with_path, const InputParameters& params);

#endif // RUNSIM_HPP