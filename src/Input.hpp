#ifndef INPUT_HPP
#define INPUT_HPP

#include "Parameters.hpp"
#include <string>
#include <iostream>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/info_parser.hpp>

std::string get_input_path(int argc, char* argv[]);
InputParameters read_parameters(const boost::property_tree::ptree& pt);

#endif // INPUT_HPP