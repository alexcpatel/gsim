#include "gsim.h"
#include "quadtree.h"

/* free all state memory */
static inline void free_state(state_t *state) {
  if (state != NULL) {
    free(state->bodies);
    free();
  }
}

/* setup state for simulation */
static int setup(state_t *state) {
  size_t num_clusters = state->num_clusters, 
         num_bodies = state->num_bodies,
         bodies_per_cluster, ci, bi;
  double cx, cy, ox, oy;
  body_t *bodies;
  partition_t *partitions;
  int thread_cnt, pi;

  /* allocate memory for bodies */
  state->bodies = calloc(num_bodies, sizeof(body_t));
  if (state->bodies == NULL) return -1;
  bodies = state->bodies;

  /* distribute random bodies by cluster */
  const double cluster_offset_scale = 
    ((double) DIST_SCALE) / ((double) num_clusters * 10);
  bodies_per_cluster = num_bodies / num_clusters;

  #pragma omp parallel for schedule(static)
  for (ci = 0; ci < num_clusters; ci++) {
    cx = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);
    cy = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);

    #pragma omp parallel for schedule(static)
    for (bi = ci * bodies_per_cluster;
         bi < (ci+1) * bodies_per_cluster
         && bi < num_bodies; bi++) {

      ox = ((double) cluster_offset_scale) *
           (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);
      oy = ((double) cluster_offset_scale) *
           (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);

      body_t b;
      b.m = ((double) INIT_MASS) +
            ((double) MASS_RANGE) *
            (((double) rand() / (double) RAND_MAX) * 2.0 - 1.0);

      b.x = cx + ox; b.vx = 0.0; b.hvx = 0.0; b.ax = 0.0;
      b.y = cy + oy; b.vy = 0.0; b.hvy = 0.0; b.ay = 0.0;
      b.work = 1;

      bodies[bi] = b;
    }
  }

  /* allocate space for partitions */
  thread_cnt = state->thread_cnt;
  state->partitions = calloc(thread_cnt, sizeof(partition_t));
  if (state->partitions == NULL) return -1;
  partitions = state->partitions;

  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    p->pbodies = calloc(num_bodies, sizeof(body_t *));
  }

  return 0;
}

/* simulate the state for one step */
static inline int step(state_t *state) {
  size_t bi, num_bodies = state->num_bodies,
         pi, thread_cnt = state->thread_cnt;
  body_t *bodies = state->bodies;
  partitions_t *partitions = state->partitions;
  quadtree_t *quadtree;
  int W, Wp;

  /* constant definitions for step loop */
  static const double epsilon2 =
    (SOFTENING_FACTOR) * (SOFTENING_FACTOR);
  static const double G = GRAV_CONSTANT;
  static const double dt = TIME_STEP;
  static const double lo_quad_bound = -1.0 * DIST_SCALE;
  static const double hi_quad_bound = 2.0  * DIST_SCALE;

  /* update positions of each body */
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);

    /* update x components */
    b->hvx = b->vx + 0.5 * b->ax * dt;
    b->x = b->x + b->hvx * dt;

    /* update y components */
    b->hvy = b->vy + 0.5 * b->ay * dt;
    b->y = b->y + b->hvy * dt;
  }

  /* build quad tree */
  quadtree = quadtree_new(lo_quad_bound, lo_quad_bound,
                          hi_quad_bound, hi_quad_bound);
  if (quadtree == NULL) return -1;

  /* insert all bodies into tree */
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    quadtree_insert(quadtree, &(bodies[bi]));
  }

  /* traverse tree to assign costs and node approximations */
  quadtree_traverse(quadtree);

  /* setup cost zone partitioning */
  W = quadtree->work;
  Wp = W / thread_cnt;
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    p->min_work = pi * Wp;
    p->max_work = pi == thread_cnt - 1 ? W : (pi + 1) * Wp;
    p->num_bodies = 0;
  }

  /* partition tree between threads by cost zone */
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    quadtree_partition(quadtree, &(partitions[pi]));
  }

  /* compute forces for each body by traversing quadtree */
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    size_t num_pbodies = p->num_pbodies;
    body_t **pbodies = p->pbodies;

    #pragma omp parallel for schedule(static)
    for (bi = 0; bi < num_pbodies; bi++) {
      body_t *b = pbodies[bi];
  
      b->work = 0;
      b->ax = 0.0;
      b->ay = 0.0;
      quadtree_aggregate_forces(quadtree, b, state->theta, epsilon2);
      b->ax *= G;
      b->ay *= G;
  
      if (b->ax < -MAX_ACCELERATION) b->ax = -MAX_ACCELERATION;
      else if (b->ax > MAX_ACCELERATION) b->ax = MAX_ACCELERATION;
      if (b->ay < -MAX_ACCELERATION) b->ay = -MAX_ACCELERATION;
      else if (b->ay > MAX_ACCELERATION) b->ay = MAX_ACCELERATION;
    }
  }
  
  /* update velocities for next step */
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);

    /* update x component */
    b->vx = b->hvx + 0.5 * b->ax * dt;

    /* update x component */
    b->vy = b->hvy + 0.5 * b->ay * dt;
  }

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
  "Usage: ./gsim <thread count> <number of clusters> <number of bodies>\n"
  "              <number of simulation steps> <theta for quadtree>\n"
  "              <random seed for initialization>\n";

int main(int argc, char **argv) {
  state_t *state;

  /* get command line arguments */
  if (argc != NUM_ARGS) {
    fprintf(stderr, usage_msg);
    return -1;
  }

  state = calloc(1, sizeof(state_t));
  if (state == NULL) return -1;

  state->thread_cnt = atoi(argv[ARG_THREAD_CNT]);
  state->num_clusters = atoi(argv[ARG_NUM_CLUSTERS]);
  state->num_bodies = atoi(argv[ARG_NUM_BODIES]);
  state->num_steps = atoi(argv[ARG_NUM_STEPS]);
  state->theta = atof(argv[ARG_THETA]);
  srand(atoi(argv[ARG_RAND_SEED]));

  fprintf(stderr, "Running with %d threads.  Max possible is %d.\n",
                  state->thread_cnt, omp_get_max_threads());

  omp_set_num_threads(state->thread_cnt);

  /* initialize simulation state */
  if (setup(state)) return -1;

  /* run simulation */
  if (simulate(state)) return -1;

  return 0;
}