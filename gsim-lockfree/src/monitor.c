#include "monitor.h"

const char* act_names[] = {
  "setup",
  "write_file",
  "update_positions",
  "build_quadtree",
  "traverse_quadtree",
  "partition_quadtree",
  "compute_forces",
  "update_velocities",
  "free_quadtree",
  "other"
};

bool   is_act[NUM_ACT];
double act_start[NUM_ACT];
double act_times[NUM_ACT];

void start_activity(act_t act) {
  if (act < 0 || act > NUM_ACT || is_act[act]) return;
  is_act[act] = true;
  act_start[act] = currentSeconds();
}

void finish_activity(act_t act) {
  if (act < 0 || act > NUM_ACT || !is_act[act]) return;
  
  act_times[act] += currentSeconds() - act_start[act];
  is_act[act] = false;
}

void print_activities(void) {
  int act;
  double total_time = 0.0;

  for (act = 0; act < NUM_ACT; act++)
    total_time += act_times[act];

  fprintf(stdout, "\n\ttotal time: %.2lf ms\n", total_time * 1000);
  fprintf(stdout, "\t--------------------------------------------------\n");

  for (act = 0; act < NUM_ACT; act++)
    fprintf(stdout, "\t%.2lf ms \t%.2lf %%    \t%s\n",
      act_times[act] * 1000,
      (act_times[act] * 100) / total_time,
      act_names[act]);

  fprintf(stdout, "\n");
}
