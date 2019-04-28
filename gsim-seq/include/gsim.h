//https://developer.download.nvidia.com/compute/cuda/1.1-Beta/x86_website/projects/nbody/doc/nbody_gems3_ch31.pdf

#ifndef _GSIM_H_
#define _GSIM_H_

#define GRAV_CONSTANT      6.67408e-11  // m^3 * kg^-1 * s^-2
#define DIST_SCALE         20.0e15      // m
#define INIT_MASS          4.0e30       // kg
#define MASS_RANGE         2.0e30       // kg
#define SOFTENING_FACTOR   1.0e10       // arbitrary small value
#define TIME_STEP          1.0e2        // s
#define MAX_ACCELERATION   5.0e6        // m/s

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

/* body struct definition */
typedef struct {
  double m;  // mass of the body

  double x;  // x coordinate of body
  double y;  // y coordinate of body

  double vx; // velocity in x direction
  double vy; // velocity in y direction

  double hvx; // half step velocity in x direction
  double hvy; // half step velocity in y direction

  double ax; // acceleration in x direction
  double ay; // acceleration in y direction
} body_t;

/* simulation state struct definition */
typedef struct {
  size_t num_clusters; // number of local clusters to divide bodies into
  size_t num_bodies; // total number of bodies
  size_t num_steps; // total number of simulation steps
  double theta; // theta for quadtree heuristic
  body_t *bodies; // array of all the data for every body
} state_t;

#endif
