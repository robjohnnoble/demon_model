#ifndef INPUT_HPP
#define INPUT_HPP

#include "parameters.hpp"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>
#include <iostream>
#include <string>

std::string getInputPath(int argc, char* argv[]);
InputParameters readParameters(const boost::property_tree::ptree& pt, const std::string& config_file_path);

#endif // INPUT_HPP