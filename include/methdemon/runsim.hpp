#ifndef RUNSIM_HPP
#define RUNSIM_HPP

#include "initialise.hpp"
#include "tumour.hpp"

#include <string>
#include <vector>

void runSim(const std::string& input_and_output_path, const std::string& config_file_with_path, const InputParameters& params);
float updateTime();

#endif // RUNSIM_HPP