//https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/nbody/doc/nbody_gems3_ch31.pdf

#ifndef _GSIM_H_
#define _GSIM_H_

#include "body.h"

#define NUM_ARGS 5

typedef enum {
  ARG_FILE_NAME = 0,
  ARG_NUM_CLUSTERS,
  ARG_NUM_BODIES,
  ARG_NUM_STEPS,
  ARG_RAND_SEED
} arg_index_t;

#define GRAV_CONSTANT 6.67408e-11

#define DIST_SCALE 20.0e15 // meters
#define INIT_MASS 2.0e30 // kg
#define MASS_RANGE 4.0e30 // kg

#define SOFTENING_FACTOR 1.0e-10
#define TIME_STEP 1.0e2 // seconds

typedef struct {
  size_t num_clusters;
  size_t num_bodies;
  size_t num_steps;
  body_t *bodies;
} state_t;

#endif
