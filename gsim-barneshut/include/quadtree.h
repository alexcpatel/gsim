#ifndef _GSIM_QUADTREE_H_
#define _GSIM_QUADTREE_H_

#include "gsim.h"

/* quadtree structure definition (192 bytes) */
typedef struct qt {
  /* open mp lock */
  omp_lock_t node_lock;

  /* tree bounding box */
  double x1; // left of bbox
  double y1; // bottom of bbox
  double x2; // right of bbox
  double y2; // top of bbox

  /* static constants to avoid computation */
  double l;  // side length of bbox
  double xm; // center of bbox
  double ym; // center of bbox

  /* accumulated tree info */
  int work;
  double m; // total mass of sub nodes
  double xc; // x center of mass
  double yc; // y center of mass

  /* leaf body */
  body_t *body;
  double bm;
  double bx;
  double by;

  /* subtrees */
  struct qt *top_left;
  struct qt *top_right;
  struct qt *bot_left;
  struct qt *bot_right;

  char pad[32];
} quadtree_t;

/* create a new quadtree */
quadtree_t *quadtree_new(double x1, double y1,
                         double x2, double y2);

/* free all of a quadtree's memory */
void quadtree_free(quadtree_t *quadtree);

/* insert a node into the quadtree */
void quadtree_insert(quadtree_t *quadtree, body_t *b);

/* post-order traversal of quadtree to compute 
 * node approximations and cumulative workloads */
void quadtree_traverse(quadtree_t *quadtree);

/* post-order traversal of the tree to get the partition
 * of quadtree nodes for the given thread */
void quadtree_partition(quadtree_t *quadtree,
                        partition_t *partition, int cur_work);

/* aggregate forces for a body from the quadtree */
void quadtree_aggregate_forces(quadtree_t *quadtree, body_t *b,
                               double theta, double epsilon2);

#endif
