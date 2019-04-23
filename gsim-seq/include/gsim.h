//https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/nbody/doc/nbody_gems3_ch31.pdf

#ifndef _GSIM_H_
#define _GSIM_H_

#include "body.h"

#define GRAV_CONSTANT      6.67408e-11  // m^3 * kg^-1 * s^-2
#define DIST_SCALE         20.0e15      // m
#define INIT_MASS          2.0e30       // kg
#define MASS_RANGE         4.0e30       // kg
#define SOFTENING_FACTOR   1.0          // arbitrary small value
#define TIME_STEP          1.0e2        // s

/* command line argument constants */
#define NUM_ARGS 6
typedef enum {
  ARG_FILE_NAME = 0,
  ARG_NUM_CLUSTERS,
  ARG_NUM_BODIES,
  ARG_NUM_STEPS,
  ARG_THETA,
  ARG_RAND_SEED
} arg_index_t;

/* simulation state struct definition */
typedef struct {
  size_t num_clusters;
  size_t num_bodies;
  size_t num_steps;
  double theta;
  body_t *bodies;
} state_t;

#endif
