#ifndef _GSIM_BODY_H_
#define _GSIM_BODY_H_

typedef struct {
  double x;  // x component of distance vector (normalized)
  double y;  // y component of distance vector (normalized)
  double mag; // magnitude of distance vector
} dist_t;

#define DIST_SCALE 20.0e15 // 20^15 meters
#define INIT_MASS 4.0e30 // 4^30 kg

typedef struct {
  double x;  // x coordinate of body
  double y;  // y coordinate of body

  double vx; // velocity in x direction
  double vy; // velocity in y direction

  double ax; // acceleration in x direction
  double ay; // acceleration in y direction

  double nax; // next acceleration in x direction
  double nay; // next acceleration in x direction

  double m;  // mass of the body
} body_t;

body_t random_body(void);

dist_t distance(body_t *b1, body_t *b2);

#endif
