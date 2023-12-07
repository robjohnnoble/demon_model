#ifndef INITIALISE_HPP
#define INITIALISE_HPP

#include "BinTrees.hpp"
#include "Parameters.hpp"
#include "Files.hpp"
#include <cmath>

DerivedParameters derived_parameters(const InputParameters& params);
float set_init_migration_rate(int K, float init_migration_rate, float A, float B, float C);
void open_output_files(const std::string& input_and_output_path, const FileStructure& layout);
void close_output_files(const std::string& input_and_output_path);

#endif // INITIALISE_HPP
