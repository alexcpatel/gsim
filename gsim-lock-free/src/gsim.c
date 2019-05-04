#include "gsim.h"
#include "quadtree.h"

/* free all state memory */
static inline void free_state(state_t *state) {
  if (state) {
    free(state->bodies);
    int pi;
    partition_t *partitions = state->partitions;
    #pragma omp parallel for schedule(static)
    for (pi = 0; pi < state->thread_cnt; pi++) {
      partition_t *p = &(partitions[pi]);
      free(p->pbodies);
    }
    free(partitions);
    free(state);
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
  if (!state->bodies) return -1;
  bodies = state->bodies;

  /* distribute random bodies by cluster */
  const double cluster_offset_scale = 
    ((double) DIST_SCALE) / ((double) num_clusters * 10);
  bodies_per_cluster = num_bodies / num_clusters + 1;

  for (ci = 0; ci < num_clusters; ci++) {
    cx = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);
    cy = ((double) DIST_SCALE) *
         ((double) rand() / (double) RAND_MAX);

    for (bi = ci * bodies_per_cluster;
         bi < (ci+1) * bodies_per_cluster; bi++) {
      if (bi < num_bodies) {
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
  }

  /* allocate space for partitions */
  thread_cnt = state->thread_cnt;
  state->partitions = calloc(thread_cnt, sizeof(partition_t));
  if (!state->partitions) return -1;
  partitions = state->partitions;

  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    p->pbodies = calloc(num_bodies, sizeof(body_t *));
  }

  return 0;
}

/* simulate the state for one step */
static inline void step(state_t *state) {
  START_ACTIVITY(OTHER);
  size_t bi, num_bodies = state->num_bodies,
         pi, thread_cnt = state->thread_cnt;
  body_t *bodies = state->bodies;
  partition_t *partitions = state->partitions;
  quadtree_t *quadtree;
  uint64_t W, Wp;

  /* constant definitions for step loop */
  static const double epsilon2 =
    (SOFTENING_FACTOR) * (SOFTENING_FACTOR);
  static const double G = GRAV_CONSTANT;
  static const double dt = TIME_STEP;
  static const double lo_quad_bound = -1.0 * DIST_SCALE;
  static const double hi_quad_bound = 2.0  * DIST_SCALE;
  FINISH_ACTIVITY(OTHER);

  /* update positions of each body */
  START_ACTIVITY(UPDATE_POSITIONS);
  //#pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);
    b->hvx = b->vx + 0.5 * b->ax * dt;
    b->x = b->x + b->hvx * dt;
    b->hvy = b->vy + 0.5 * b->ay * dt;
    b->y = b->y + b->hvy * dt;
  }
  FINISH_ACTIVITY(UPDATE_POSITIONS);

  /* initialize new quadtree */
  START_ACTIVITY(BUILD_QUADTREE);
  quadtree = quadtree_new(lo_quad_bound, lo_quad_bound,
                          hi_quad_bound, hi_quad_bound);

  /* insert all bodies into tree */
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    quadtree_insert(quadtree, &(bodies[bi]));
  }
  FINISH_ACTIVITY(BUILD_QUADTREE);

  /* traverse tree to assign costs and node approximations */
  START_ACTIVITY(TRAVERSE_QUADTREE);
  quadtree_traverse(quadtree);
  FINISH_ACTIVITY(TRAVERSE_QUADTREE);

  /* setup cost zone partitioning */
  START_ACTIVITY(PARTITION_QUADTREE);
  W = quadtree->work;
  Wp = W / thread_cnt;
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    p->min_work = pi * Wp;
    p->max_work = pi == thread_cnt - 1 ? W : (pi + 1) * Wp - 1;
    if (p->max_work < p->min_work) p->max_work = p->min_work;
    p->num_pbodies = 0;
  }

  /* partition tree between threads by cost zone */
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    quadtree_partition(quadtree, &(partitions[pi]), 0);
  }
  FINISH_ACTIVITY(PARTITION_QUADTREE);

  /* compute forces for each body by traversing quadtree */
  START_ACTIVITY(COMPUTE_FORCES);
  #pragma omp parallel for schedule(static)
  for (pi = 0; pi < thread_cnt; pi++) {
    partition_t *p = &(partitions[pi]);
    size_t num_pbodies = p->num_pbodies;
    body_t **pbodies = p->pbodies;

    size_t pbi;
    for (pbi = 0; pbi < num_pbodies; pbi++) {
      body_t *b = pbodies[pbi];
  
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
  FINISH_ACTIVITY(COMPUTE_FORCES);
  
  /* update velocities for next step */
  START_ACTIVITY(UPDATE_VELOCITIES);
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);
    b->vx = b->hvx + 0.5 * b->ax * dt;
    b->vy = b->hvy + 0.5 * b->ay * dt;
  }
  FINISH_ACTIVITY(UPDATE_VELOCITIES);

  START_ACTIVITY(FREE_QUADTREE);
  quadtree_free(quadtree);
  FINISH_ACTIVITY(FREE_QUADTREE);
}

/* writes the header line of the simulation state to stdout */
static inline void write_header(state_t *state) {
  START_ACTIVITY(WRITE_FILE);
  fprintf(state->output_file, "%zu %zu\n",
          state->num_bodies, state->num_steps);
  FINISH_ACTIVITY(WRITE_FILE);
}

/* writes the current simulation state to stdout */
static inline void write_state(state_t *state, size_t si) {
  START_ACTIVITY(WRITE_FILE);
  size_t bi;
  body_t *b;

  for (bi = 0; bi < state->num_bodies; bi++) {
    b = &(state->bodies[bi]);
    fprintf(state->output_file, "%zu %zu %f %f %f\n",
            si, bi, b->m, b->x, b->y);
  }
  FINISH_ACTIVITY(WRITE_FILE);
}

/* simulate the state for a number of steps */
static void simulate(state_t *state) {
  size_t si;
  bool has_output_file = !!state->output_file;

  if (has_output_file) write_header(state);
  for (si = 0; si < state->num_steps; si++) {
    if (has_output_file) write_state(state, si);
    step(state);
  }
  if (has_output_file) write_state(state, state->num_steps);
}

static const char *usage_msg =
  "Usage: ./gsim <thread count> <number of clusters> <number of bodies>\n"
  "              <number of simulation steps> <theta for quadtree>\n"
  "              <random seed for initialization> [output file]\n";

int main(int argc, char **argv) {
  state_t *state;
  START_ACTIVITY(OTHER);

  /* get command line arguments */
  if (argc < MIN_ARGS || argc > MAX_ARGS) {
    fprintf(stderr, usage_msg);
    return -1;
  }

  if (!(state = calloc(1, sizeof(state_t)))) return -1;
  state->thread_cnt = atoi(argv[ARG_THREAD_CNT]);
  state->num_clusters = atoi(argv[ARG_NUM_CLUSTERS]);
  state->num_bodies = atoi(argv[ARG_NUM_BODIES]);
  state->num_steps = atoi(argv[ARG_NUM_STEPS]);
  state->theta = atof(argv[ARG_THETA]);
  srand(atoi(argv[ARG_RAND_SEED]));
  state->output_file = argc == MAX_ARGS ? 
    fopen(argv[ARG_OUT_FILE], "w+") : NULL;

  fprintf(stdout, "Running with %d %s. Max possible is %d.\n",
                  state->thread_cnt,
                  state->thread_cnt == 1 ? "thread" : "threads",
                  omp_get_max_threads());

  omp_set_num_threads(state->thread_cnt);
  FINISH_ACTIVITY(OTHER);

  // quadtree_t *test = quadtree_new(0.0,0.0,0.0,0.0);
  // test->is_leaf = 0x7;
  // test->has_body = 0x2;
  // fprintf(stdout, "%u, %u\n", (unsigned int) test->is_leaf, (unsigned int) test->has_body);
  // fprintf(stdout, "%x\n", (unsigned int) *((uint32_t *) &(test->is_leaf)));
  // quadtree_free(test);
  // exit(0);

  /* initialize simulation state */
  START_ACTIVITY(SETUP);
  if (setup(state)) return -1;
  FINISH_ACTIVITY(SETUP);

  /* run simulation */
  simulate(state);

  /* free and close resources */
  START_ACTIVITY(WRITE_FILE);
  if (state->output_file) fclose(state->output_file);
  FINISH_ACTIVITY(WRITE_FILE);

  START_ACTIVITY(OTHER);
  free_state(state);
  FINISH_ACTIVITY(OTHER);

  /* print monitor information */
  PRINT_ACTIVITIES();
  return 0;
}