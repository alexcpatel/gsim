#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "body.h"
#include "gsim.h"

/* setup state for simulation */
static int setup(state_t *state) {
  size_t num_clusters = state->num_clusters, 
         num_bodies = state->num_bodies,
         bodies_per_cluster, ci, bi;
  double cx, cy, ox, oy;
  body_t *bodies, b;

  /* allocate memory for bodies */
  state->bodies = calloc(num_bodies, sizeof(body_t));
  if (state->bodies == NULL) return -1;
  bodies = state->bodies;

  /* distribute random bodies by cluster */
  const double cluster_offset_scale = 
    ((double) DIST_SCALE) / ((double) num_clusters * 10);
  bodies_per_cluster = num_bodies / num_clusters;

  for (ci = 0; ci < num_clusters; ci++) {
    cx = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);
    cy = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);

    for (bi = ci * bodies_per_cluster;
         bi < (ci+1) * bodies_per_cluster
         && bi < num_bodies; bi++) {

      ox = ((double) cluster_offset_scale) *
           (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);
      oy = ((double) cluster_offset_scale) *
           (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);

      b.x = cx + ox;
      b.y = cy + oy;
    
      b.vx = 0.0; b.ax = 0.0; b.nax = 0.0;
      b.vy = 0.0; b.ay = 0.0; b.nay = 0.0;

      b.m = ((double) INIT_MASS) +
            ((double) MASS_RANGE) *
            (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);

      bodies[bi] = b;


    }
  }

  return 0;
}

/* simulate the state for one step */
static inline int step(state_t *state) {
  size_t b1i, b2i, num_bodies = state->num_bodies;
  double ax, ay, mult, denom, denom_sq;
  body_t *bodies = state->bodies;
  body_t *b1, *b2;
  dist_t dist;

  const double epsilon_sq = (SOFTENING_FACTOR) * (SOFTENING_FACTOR);
  const double G = GRAV_CONSTANT;
  const double dt = TIME_STEP;
  const double dt2 = (dt*dt) / 2.0;

  /* compute forces */
  for (b1i = 0; b1i < num_bodies; b1i++) {
    b1 = &(bodies[b1i]);
    
    ax = 0.0;
    ay = 0.0;

    /* add contributions from all other bodies */
    for (b2i = 0; b2i < num_bodies; b2i++) {
      if (b1i == b2i) continue;

      b2 = &(bodies[b2i]);
      dist = distance(b1, b2);

      denom = dist.mag*dist.mag + epsilon_sq;
      denom_sq = denom*denom;
      denom = (denom * denom_sq) / denom_sq;
      mult = b2->m / denom;

      ax += mult * dist.x;
      ay += mult * dist.y;
    }

    b1->nax = G * ax;
    b1->nay = G * ay;
  }

  //fprintf(stderr, "FINISHED: computing forces!\n");

  /* update positions */
  for (b1i = 0; b1i < num_bodies; b1i++) {
    b1 = &(bodies[b1i]);

    /* update x components */
    b1->x = b1->x + b1->vx * dt + b1->ax * dt2;
    b1->vx = b1->vx + ((b1->nax + b1->ax) / 2.0) * dt;
    b1->ax = b1->nax;

    /* update y components */
    b1->y = b1->y + b1->vy * dt + b1->ay * dt2;
    b1->vy = b1->vy + ((b1->nay + b1->ay) / 2.0) * dt;
    b1->ay = b1->nay;
  }

  //fprintf(stderr, "FINISHED: updating positions!\n");

  return 0;
}

/* writes the header line of the simulation state to stdout */
static inline void write_header(state_t *state) {
  fprintf(stdout, "%zu %zu\n", state->num_bodies, state->num_steps);
}

/* writes the current simulation state to stdout */
static inline void write_state(state_t *state, size_t si) {
  size_t bi;
  body_t *b;

  for (bi = 0; bi < state->num_bodies; bi++) {
    b = &(state->bodies[bi]);
    fprintf(stdout, "%zu %zu %f %f %f\n",
            si, bi, b->m, b->x, b->y);
  }
}

/* simulate the state for a number of steps */
static int simulate(state_t *state) {
  size_t si;

  write_header(state);

  for (si = 0; si < state->num_steps; si++) {
    write_state(state, si);
    if (step(state)) return -1;
  }

  write_state(state, state->num_steps);

  return 0;
}

static const char *usage_msg =
  "Usage: ./gsim <number of clusters> <number of bodies> "
  "<number of steps> <random seed for initialization>\n";

int main(int argc, char **argv) {
  state_t *state;

  /* get command line arguments */
  if (argc != NUM_ARGS) {
    fprintf(stderr, usage_msg);
    return -1;
  }

  state = calloc(1, sizeof(state_t));
  if (state == NULL) return -1;

  state->num_clusters = atoi(argv[ARG_NUM_CLUSTERS]);
  state->num_bodies = atoi(argv[ARG_NUM_BODIES]);
  state->num_steps = atoi(argv[ARG_NUM_STEPS]);
  srand(atoi(argv[ARG_RAND_SEED]));

  /* initialize simulation state */
  if (setup(state)) return -1;

  fprintf(stderr, "finished setting up!\n");

  /* run simulation */
  if (simulate(state)) return -1;

  fprintf(stderr, "finished simulation!\n");

  return 0;
}