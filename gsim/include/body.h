#ifndef _GSIM_BODY_H_
#define _GSIM_BODY_H_

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

#endif
