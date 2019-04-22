//https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/nbody/doc/nbody_gems3_ch31.pdf

#ifndef _GSIM_H_
#define _GSIM_H_

#include "body.h"

#define NUM_ARGS 3

typedef enum {
  ARG_FILE_NAME = 0,
  ARG_NUM_BODIES,
  ARG_NUM_STEPS
} arg_index_t;

#define SOFTENING_FACTOR 1.0e-10
#define GRAV_CONSTANT 6.67408e-11
#define TIME_STEP 1.0e2

#define DIST_THRESH 20.0

typedef struct {
  size_t num_bodies;
  size_t num_steps;
  body_t *bodies;
} state_t;

#endif
