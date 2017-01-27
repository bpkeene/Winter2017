#!/bin/bash

rm -f *aux *log *gz *.out *exe
g++ simulator.cpp random_mars.cpp -lfftw3 -lm -o MC_exe;
