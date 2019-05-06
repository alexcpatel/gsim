# gsim
Parallel Galaxy Simulation with the Barnes-Hut Algorithm using OpenMP

## Overview

* The `gsim-bad` directory contains an optimized parallel implementation of a naive all-pairs O(n^2) algorithm to act as a baseline. 
* The `gsim-barneshut` directory contains an optimized parallel implementation of the Barnes-Hut Algorithm with a fine-grained locking quadtree for concurrent inserts.
* The `gsim-lockfree` directory contains an optimized parallel implementation of the Barnes-Hut Algorithm with a lock-free quadtree for concurrent inserts.

The above project directories have a Makefile (run `make` in the project directory) to build the respective projects. Run the generated binary without arguments to view usage. Optionally, run the executable with a specified output file to generate data for the visualizer.

The `gviz` directory contains a program to accept an output file from either of the two versions of the implementation and display a visualization. This directory also contains a Makefile, and the visualizer requires `glfw3` installed on your system.

## Benchmarking

The benchmarking script `benchmark.py` runs all of our implementations on varying parameters and writes the results into an output file. The benchmark also writes trace results for each run performed.

We benchmark over `gsim-bad`, `gsim-barneshut`, and `gsim-lockfree`, with theta in `[0.1, 0.3, 0.5]`. For each of these cases we have 3 benchmarks:

 * Benchmark *A*: 1-to-1, an equal number of clusters and bodies in the simulation
 * Benchmark *B*: sqrt, the number of clusters is the square root of the number of bodies
 * Benchmark *C*: single, all bodies are initialized in a single cluster

The benchmark script supports the following optional flags:

 * `-t`: Maximum number of threads to perform runs up to (machine dependent) (default: `16`)
 * `-b`: Total number of simulation bodies for each run (default: `1000`)
 * `-s`: Total number of simulation steps for each run (default: `100`)
 * `-r`: Number of runs we take the best result from for each benchmark
 * `-f`: Path to directory traces should be written to (default: `./traces`)
 * `-o`: Output file of benchmark results (default: `./benchmark.txt`)
