#!/opt/homebrew/zsh

g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ main.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Clone.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Deme.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Genotype.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Initialise.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Input.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Output.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ RunSim.cpp
g++ -g -c -I/opt/homebrew/Cellar/boost/1.82.0_1/include/ Tumour.cpp

g++ -v -o methdemon main.o Clone.o Deme.o Genotype.o Initialise.o Input.o RunSim.o Tumour.o -lm
