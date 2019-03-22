# demon
Demon (deme-based oncology model) is a flexible framework for modelling intra-tumour population genetics with varied spatial structures and modes of cell dispersal.

## Prerequisites

The program is written in C++ (mostly plain C); it requires gnuplot and the Boost Property Tree C++ library.

## Running a simulation

Copy the following three files into a folder:

* `demon.cpp` (model code);
* `demon.h` (header file);
* `init_conf_file.dat` (configuration file);
* `compile_execute.sh` (shell script to compile and execute the model).

In a terminal window, navigate to your folder and type `./compile_execute.sh`.

## Configuration file

The configuration file sets the following model parameters:

Spatial structure
* `int log2_deme_carrying_capacity`: log2 of deme carrying capacity

Dispersal
* `float migration_type`: type of cell dispersal (0 = invasion; 2 = deme fission)
* `float init_migration_rate`: initial invasion rate or deme fission rate
* `int migration_edge_only`: whether dispersal is limited to the tumour boundary (0 = no; 1 = yes)
* `int migration_rate_scales_with_K`: whether to divide dispersal rates by the square root of the deme carrying capacity (0 = no; 1= yes)

Mutation rates
* `float mu_driver_birth`: driver mutation rate per cell division (for drivers affecting division rate)
* `float mu_passenger`: passenger mutation rate per cell division
* `float mu_driver_migration`: driver mutation rate per cell division (for drivers affecting dispersal rate)
* `int passenger_pop_threshold`: population size at which passenger mutations stop occurring (-1 = no threshold)

Fitness effects
* `float normal_birth_rate`: normal cell division rate, relative to tumour cell division rate (-1 = no normal cells)
* `float baseline_death_rate`: baseline death rate, independent of deme population size
* `float s_passenger`: deleterious effect on division rate per passenger mutation
* `float s_driver_birth`: beneficial effect on division rate per driver mutation (for drivers affecting division rate)
* `float s_driver_migration`: beneficial effect on dispersal rate per driver mutation (for drivers affecting dispersal rate)
* `float max_relative_birth_rate`: maximum division rate, relative to initial division rate
* `float max_relative_migration_rate`: maximum dispersal rate, relative to initial dispersal rate

Non-biological parameters
* `int seed`: seed for pseudo-random number generator
* `int init_pop`: initial tumour cell population
* `int max_time`: max elapsed time before program halts, in seconds
* `int max_pop`: max tumour cell population before program halts
* `int max_generations`: max cell generations before program halts
* `int matrix_max`: max number of genotypes before program halts (larger value => more allocated memory)
* `int write_grid`: whether to generate plots of grid states during execution (can significantly increase running time)
* `int write_clones_file`: whether to write contents of all clones (can be a very large file; typically unnecessary)
* `int write_demes_file`: whether to write contents of all demes (can be a very large file; typically unnecessary)
* `int record_matrix`: whether to record genetic distance matrix for all genotypes (can be a very large file; typically unnecessary)
* `int write_phylo`: whether to write phylogeny for all genotypes (can be a very large file; typically unnecessary)
* `int calculate_total_diversity`: whether to calculate diversity across all genotypes (can be computationally expensive)
* `int biopsy_size_per_sample`: max number of cells per biopsy sample (reserved for future applications)

