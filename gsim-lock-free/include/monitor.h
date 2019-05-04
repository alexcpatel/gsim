#ifndef ACTIVITY_MONITOR_H
#define ACTIVITY_MONITOR_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <omp.h>

#include "cycletimer.h"

/* monitor activity macros */
#define DEBUG 1

#if     DEBUG
#define START_ACTIVITY(x)     do { start_activity(x);  } while(0)
#define FINISH_ACTIVITY(x)    do { finish_activity(x); } while(0)
#define PRINT_ACTIVITIES(x)   do { print_activities(); } while(0)
#else
#define START_ACTIVITY(x)
#define FINISH_ACTIVITY(x)
#define PRINT_ACTIVITIES(x)
#endif

/* activity definitions */
typedef enum {
  SETUP = 0,
  WRITE_FILE,
  UPDATE_POSITIONS,
  BUILD_QUADTREE,
  TRAVERSE_QUADTREE,
  PARTITION_QUADTREE,
  COMPUTE_FORCES,
  UPDATE_VELOCITIES,
  FREE_QUADTREE,
  OTHER,
  NUM_ACT
} act_t;

/* monitor state structures */
extern const char* act_names[];
extern bool   is_act[NUM_ACT];
extern double act_start[NUM_ACT];
extern double act_times[NUM_ACT];

/* monitor interface */
void start_activity(act_t act);
void finish_activity(act_t act);
void print_activities(void);

#endif
