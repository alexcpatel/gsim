#ifndef _GSIM_H_
#define _GSIM_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>

#include "monitor.h"

/* general constants for gsim */
#define GRAV_CONSTANT      6.67408e-11  // m^3 * kg^-1 * s^-2
#define DIST_SCALE         20.0e15      // m
#define INIT_MASS          4.0e30       // kg
#define MASS_RANGE         2.0e30       // kg
#define SOFTENING_FACTOR   1.0          // arbitrary small value
#define TIME_STEP          1.0e2        // s
#define MAX_ACCELERATION   7.0e6        // m * s^-2

/* command line argument constants */
#define MIN_ARGS 6
#define MAX_ARGS 7
typedef enum {
  ARG_FILE_NAME = 0,
  ARG_THREAD_CNT,
  ARG_NUM_CLUSTERS,
  ARG_NUM_BODIES,
  ARG_NUM_STEPS,
  ARG_RAND_SEED,
  ARG_OUT_FILE
} arg_index_t;

/* body struct definition (128 bytes) */
typedef struct bt {
  double m;  // mass of the body

  double x;  // x coordinate of body
  double y;  // y coordinate of body

  double vx; // velocity in x direction
  double vy; // velocity in y direction

  double hvx; // half step velocity in x direction
  double hvy; // half step velocity in y direction

  double ax; // acceleration in x direction
  double ay; // acceleration in y direction

  char pad[56];
} body_t;

/* simulation state struct definition (64 bytes) */
typedef struct {
  FILE *output_file; // output file for visualizer
  body_t *bodies; // array of all the data for every body
  size_t num_clusters; // number of local clusters to divide bodies into
  size_t num_bodies; // total number of bodies
  size_t num_steps; // total number of simulation steps  
  int thread_cnt; // number of openmp threads
  char pad[20];
} state_t;

#endif
