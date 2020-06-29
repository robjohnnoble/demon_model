#!/bin/bash

mkdir -p BirthRateGrids
mkdir -p DriverGrids
mkdir -p PopsGrids
mkdir -p PassengersGrids
mkdir -p MigrationRateGrids

if g++ -o demon demon.cpp -I/usr/local/include/boost/ -lm; then #compile
	echo "Compiled"
	./demon  ./ init_conf_file.dat #execute
	printf "\nEnd of execution\n"
else
	echo "Failed to compile"
fi