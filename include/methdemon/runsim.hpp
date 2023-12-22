#ifndef RUNSIM_HPP
#define RUNSIM_HPP

#include "initialise.hpp"
#include "output.hpp"
#include "tumour.hpp"

#include <chrono>
#include <string>
#include <vector>

void runSim(const std::string& input_and_output_path, const std::string& config_file_with_path, const InputParameters& params);
float calculateTime(Tumour& tumour);

#endif // RUNSIM_HPP