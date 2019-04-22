# gsim-mpi
Parallel Galaxy Simulation with the Barnes Hut Algorithm using MPI

## Overview

The `gsim` directory contains a sequential version of the algorithm.
The `gsim-mpi` directory contains a parallel version of the algorithm.

Both of the above project directories have a makefile to build the respective projects. Run the executable and pipe stdout to a text file to generate the output.

The `gviz` directory contains a program to accept an output file from either of the two versions of the implementation and display a visualization.
