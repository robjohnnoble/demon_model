#ifndef INITIALISE_HPP
#define INITIALISE_HPP

//#include "BinTrees.hpp"
#include "parameters.hpp"
#include "macros.hpp"
#include "tumour.hpp"
#include <cmath>
#include <string>

DerivedParameters deriveParameters(const InputParameters& params);
int calculateMaxGenerations(float fission_rate, int deme_carrying_capacity);
// TO DO: files for output - most functions needed

#endif // INITIALISE_HPP
