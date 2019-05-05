#include "gsim.h"

/* free all state memory */
static inline void free_state(state_t *state) {
  if (state) {
    free(state->bodies);
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

  /* allocate memory for bodies */
  state->bodies = calloc(num_bodies, sizeof(body_t));
  if (state->bodies == NULL) return -1;
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
  
        bodies[bi] = b;
      }
    }
  }

  return 0;
}

/* simulate the state for one step */
static inline void step(state_t *state) {
  START_ACTIVITY(OTHER);
  size_t bi, num_bodies = state->num_bodies;
  body_t *bodies = state->bodies;

  /* constant definitions for step loop */
  static const double epsilon2 =
    (SOFTENING_FACTOR) * (SOFTENING_FACTOR);
  static const double G = GRAV_CONSTANT;
  static const double dt = TIME_STEP;
  FINISH_ACTIVITY(OTHER);

  /* update positions of each body */
  START_ACTIVITY(UPDATE_POSITIONS);
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);
    b->hvx = b->vx + 0.5 * b->ax * dt;
    b->x = b->x + b->hvx * dt;
    b->hvy = b->vy + 0.5 * b->ay * dt;
    b->y = b->y + b->hvy * dt;
  }
  FINISH_ACTIVITY(UPDATE_POSITIONS);

  /* compute forces for each body by traversing quadtree */
  START_ACTIVITY(COMPUTE_FORCES);
  #pragma omp parallel for schedule(static)
  for (bi = 0; bi < num_bodies; bi++) {
    body_t *b = &(bodies[bi]);

    b->ax = 0.0;
    b->ay = 0.0;

    double x = b->x;
    double y = b->y;

    size_t obi;
    for (obi = 0; obi < num_bodies; obi++) {
      if (obi == bi) continue;

      body_t *ob = &(bodies[obi]);
      double dx = ob->x - x;
      double dy = ob->y - y;
      double dmag = sqrt(dx*dx + dy*dy); 
      double mult = ob->m / (dmag * dmag + epsilon2);
      
      b->ax += mult * dx;
      b->ay += mult * dy;
    }

    b->ax *= G;
    b->ay *= G;
  
    if (b->ax < -MAX_ACCELERATION) b->ax = -MAX_ACCELERATION;
    else if (b->ax > MAX_ACCELERATION) b->ax = MAX_ACCELERATION;
    if (b->ay < -MAX_ACCELERATION) b->ay = -MAX_ACCELERATION;
    else if (b->ay > MAX_ACCELERATION) b->ay = MAX_ACCELERATION;
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
  bool has_output_file = state->output_file != NULL;

  if (has_output_file) write_header(state);
  for (si = 0; si < state->num_steps; si++) {
    if (has_output_file) write_state(state, si);
    step(state);
  }
  if (has_output_file) write_state(state, state->num_steps);
}

static const char *usage_msg =
  "Usage: ./gsim-bad <thread count> <number of clusters> <number of bodies>\n"
  "                  <number of simulation steps> <random seed for initialization>\n"
  "                  [output file]\n";

int main(int argc, char **argv) {
  state_t *state;
  START_ACTIVITY(OTHER);

  /* get command line arguments */
  if (argc < MIN_ARGS || argc > MAX_ARGS) {
    fprintf(stderr, usage_msg);
    return -1;
  }

  state = calloc(1, sizeof(state_t));
  if (state == NULL) return -1;

  state->thread_cnt = atoi(argv[ARG_THREAD_CNT]);
  state->num_clusters = atoi(argv[ARG_NUM_CLUSTERS]);
  state->num_bodies = atoi(argv[ARG_NUM_BODIES]);
  state->num_steps = atoi(argv[ARG_NUM_STEPS]);
  srand(atoi(argv[ARG_RAND_SEED]));
  state->output_file = argc == MAX_ARGS ? 
    fopen(argv[ARG_OUT_FILE], "w+") : NULL;

  fprintf(stderr, "Running with %d %s. Max possible is %d.\n",
                  state->thread_cnt,
                  state->thread_cnt == 1 ? "thread" : "threads",
                  omp_get_max_threads());

  omp_set_num_threads(state->thread_cnt);
  FINISH_ACTIVITY(OTHER);

  /* initialize simulation state */
  START_ACTIVITY(SETUP);
  if (setup(state)) return -1;
  FINISH_ACTIVITY(SETUP);

  /* run simulation */
  simulate(state);

  /* free and close resources */
  START_ACTIVITY(WRITE_FILE);
  if (state->output_file != NULL) fclose(state->output_file);
  FINISH_ACTIVITY(WRITE_FILE);

  START_ACTIVITY(OTHER);
  free_state(state);
  FINISH_ACTIVITY(OTHER);

  /* print monitor information */
  PRINT_ACTIVITIES();
  return 0;
}