#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "body.h"

/* create a random body */
body_t random_body(void) {
  body_t b;

  b.x = (double) DIST_SCALE *
        ((double) rand() / (double) RAND_MAX);
  b.y = (double) DIST_SCALE *
        ((double) rand() / (double) RAND_MAX);

  b.vx = 0.0;
  b.vy = 0.0;

  b.ax = 0.0;
  b.ay = 0.0;

  b.nax = 0.0;
  b.nay = 0.0;

  b.m = INIT_MASS;

  return b;
}

/* compute */
dist_t distance(body_t *b1, body_t *b2) {
  double dx = b2->x - b1->x;
  double dy = b2->y - b1->y;
  double mag = sqrt(dx*dx + dy*dy);

  dist_t dist;
  dist.x = dx;
  dist.y = dy;
  dist.mag = mag;

  return dist;
}