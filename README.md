# demon
Demon (deme-based oncology model) is a flexible framework for modelling intra-tumour population genetics with varied spatial structures and modes of cell dispersal.

## Prerequisites

The program requires gnuplot and the Boost Property Tree C++ library and a C++ compiler (such as GCC).

## Running a simulation

Copy the following three files into a folder:

* `demon.cpp` (model code);
* `demon.h` (header file);
* `init_conf_file.dat` (configuration file);
* `compile_execute.sh` (shell script to compile and execute the model).

In a terminal window, navigate to your folder and type `./compile_execute.sh`.
