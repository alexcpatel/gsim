#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "body.h"

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