# gsim
Parallel Galaxy Simulation with the Barnes-Hut Algorithm using OpenMP

## Overview

The `gsim-bad` directory contains an implementation of a naive O(n^2) algorithm. 
The `gsim-seq` directory contains a baseline sequential implementation of the algorithm.
The `gsim-omp` directory contains a parallel implementation of the algorithm.

The above project directories have a Makefile (run `make` in the project directory) to build the respective projects. Run the executable with a specified output file to generate data for the visualizer.

The `gviz` directory contains a program to accept an output file from either of the two versions of the implementation and display a visualization. This directory also contains a Makefile, and the visualizer requires `glfw3` installed on your system.

## Benchmarking
The benchmarking script `benchmark.py` runs our OpenMP implementation over multiple parameters and compares the outputs. Parameters that are varied between benchmarks include: number of threads, total number of bodies, and number of local clusters bodies are divided into. Other parameters which can be specified include: random seed for initialization, number of simulation steps, output file to view visualizer trace. Each benchmark takes the maximum out of some number of runs (default 10), so we can get an optimistic estimate of the performance for that benchmark.

Run `./benchmark.py -h` to examine optional flags. The script will run with default settings if no flags are specified.
