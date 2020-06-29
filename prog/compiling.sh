#!/bin/bash

#mkdir -p BirthRateGrids
#mkdir -p DriverGrids
#mkdir -p PopsGrids
#mkdir -p PassengersGrids

g++ -o demon demon.cpp -I/usr/local/include/boost/ -lm #compile

printf "\nCode has compiled\n"