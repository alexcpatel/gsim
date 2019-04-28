#ifndef _GSIM_QUADTREE_H_
#define _GSIM_QUADTREE_H_

#include "gsim.h"

typedef struct qt {
  /* tree bounding box */
  double x1; // left of bbox
  double y1; // bottom of bbox
  double x2; // right of bbox
  double y2; // top of bbox

  /* static constants to avoid computation */
  double l;  // side length of bbox
  double xm; // center of bbox
  double ym; // center of bbox

  /* tree info */
  double m;  // total mass of sub nodes
  double xc; // x center of mass
  double yc; // y center of mass

  /* leaf node data */
  bool has_node; // does this tree have data?
  double nm; // mass of node
  double nx; // x of node
  double ny; // y of node

  /* subtrees */
  struct qt *top_left;
  struct qt *top_right;
  struct qt *bot_left;
  struct qt *bot_right;
} quadtree_t;

/* create a new quadtree */
quadtree_t *quadtree_new(double x1, double y1,
                         double x2, double y2);

/* insert a node into the quadtree */
int quadtree_insert(quadtree_t *quadtree, double m,
                    double x, double y);

/* aggregate forces for a body from the quadtree */
int aggregate_forces(quadtree_t *quadtree, body_t *b,
                     double theta, double epsilon2);

#endif
