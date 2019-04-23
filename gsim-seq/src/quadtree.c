#include <stdlib.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

#include "quadtree.h"

/* create a new quadtree */
quadtree_t *quadtree_new(double x1, double y1, 
                         double x2, double y2) {
  quadtree_t *quadtree = calloc(1, sizeof(quadtree_t));
  if (quadtree == NULL) return NULL;

  quadtree->x1 = x1;
  quadtree->y1 = y1;
  quadtree->x2 = x2;
  quadtree->y2 = y2;

  quadtree->l  = x2 - x1;
  quadtree->xm = 0.5 * (x1 + x2);
  quadtree->ym = 0.5 * (y1 + y2);

  quadtree->m  = 0.0;
  quadtree->xc = quadtree->xm;
  quadtree->yc = quadtree->ym;

  quadtree->has_node = false;
  quadtree->nm = 0.0;
  quadtree->nx = 0.0;
  quadtree->ny = 0.0;

  quadtree->top_left  = NULL;
  quadtree->top_right = NULL;
  quadtree->bot_left  = NULL;
  quadtree->bot_right = NULL;

  return quadtree;
}

/* check if a position is in the bounding box of the quadtree*/
static inline bool in_boundary(quadtree_t *quadtree,
                               double x, double y) {
  return quadtree != NULL &&
         x >= quadtree->x1 &&
         x <= quadtree->x2 &&
         y >= quadtree->y1 &&
         y <= quadtree->y2;
}

/* check if the quadtree is a leaf (has no subtrees) */
static inline bool is_leaf(quadtree_t *quadtree) {
  return quadtree != NULL &&
         quadtree->bot_left == NULL &&
         quadtree->top_left == NULL &&
         quadtree->bot_right == NULL &&
         quadtree->top_right == NULL;
}

/* get the next subtree for the position */
static inline quadtree_t *get_subtree(quadtree_t *quadtree,
                                      double x, double y) {

  /* choose the next subtree */
  if (x < quadtree->xm) {
    if (y < quadtree->ym) return quadtree->bot_left;
    else return quadtree->top_left;
  } else {
    if (y < quadtree->ym) return quadtree->bot_right;
    else return quadtree->top_right;
  }

  return NULL;
}

/* insert a body into the quadtree */
int quadtree_insert(quadtree_t *quadtree, double m,
                     double x, double y) {
  double x1, y1, x2, y2, xm, ym, nm, nx, ny;
  quadtree_t *q1, *q2;

  /* check if body is in the quadtree */
  if (!in_boundary(quadtree, x, y)) return 0;

  x1 = quadtree->x1;
  y1 = quadtree->y1;
  x2 = quadtree->x2;
  y2 = quadtree->y2;
  xm = quadtree->xm;
  ym = quadtree->ym;

  /* update tree mass and center of mass */
  quadtree->xc = quadtree->xc * quadtree->m + x * m;
  quadtree->yc = quadtree->yc * quadtree->m + y * m;
  quadtree->m += m;
  quadtree->xc /= quadtree->m;
  quadtree->yc /= quadtree->m;

  /* attempt to insert body at leaf */
  if (is_leaf(quadtree)) {

    /* immediately insert if quadtree has no tree */
    if (!quadtree->has_node) {
      quadtree->has_node = true;
      quadtree->nm = m;
      quadtree->nx = x;
      quadtree->ny = y;
      return 0;
    } else {

      /* there is a collision, must allocate subtrees*/
      quadtree->bot_left = quadtree_new(x1, y1, xm, ym);
      quadtree->top_left = quadtree_new(x1, ym, xm, y2);
      quadtree->bot_right = quadtree_new(xm, y1, x2, ym);
      quadtree->top_right = quadtree_new(xm, ym, x2, y2);

      if (quadtree->bot_left  == NULL ||
          quadtree->top_left  == NULL ||
          quadtree->bot_right == NULL ||
          quadtree->top_right == NULL) return -1;

      /* remove node and re-insert */
      quadtree->has_node = false;
      nm = quadtree->nm;
      nx = quadtree->nx;
      ny = quadtree->nx;

      q1 = get_subtree(quadtree, nx, ny);
      if (q1 == NULL) return -1;
      if (quadtree_insert(q1, nm, nx, ny)) return -1;
    }
  }

  /* recurse to subtrees */
  q2 = get_subtree(quadtree, x, y);
  if (q2 == NULL) return -1;
  if (quadtree_insert(q2, m, x, y)) return -1;

  return 0;
}

/* aggregate forces for a body using the tree */
int aggregate_forces(quadtree_t *quadtree, body_t *b,
                     double theta, double epsilon2) {
  double x = b->x, y = b->y, dx, dy, dmag,
         denom, denom2, mult;
  int err;

  /* check if the condition is satisfied */
  dx = quadtree->xc - x;
  dy = quadtree->yc - y;
  dmag = sqrt(dx*dx + dy*dy);

  if (quadtree->l/dmag < theta) {
    fprintf(stderr, "segfault\n");

    /* add tree contribution */
    denom = dmag*dmag + epsilon2;
    denom2 = denom*denom;
    denom = (denom * denom2) / denom2;
    mult = quadtree->m / denom;

    b->ax += mult * dx;
    b->ay += mult * dy;

    return 0;
  }

  /* add node contribution if it exists */
  if (quadtree->has_node) {
    dx = quadtree->nx - x;
    dy = quadtree->ny - y;
    dmag = sqrt(dx*dx + dy*dy);

    denom = dmag*dmag + epsilon2;
    denom2 = denom*denom;
    denom = (denom * denom2) / denom2;
    mult = quadtree->nm / denom;

    b->ax += mult * dx;
    b->ay += mult * dy;

    return 0;
  }

  /* recurse to aggregate subtrees */
  err = 0;

  if (quadtree->bot_left != NULL) {
    err |= aggregate_forces(quadtree->bot_left, b, 
                            theta, epsilon2);
  }

  if (quadtree->top_left != NULL) {
    err |= aggregate_forces(quadtree->top_left, b, 
                            theta, epsilon2);
  }

  if (quadtree->bot_right != NULL) {
    err |= aggregate_forces(quadtree->bot_right, b, 
                            theta, epsilon2);
  }

  if (quadtree->top_right != NULL) {
    err |= aggregate_forces(quadtree->top_right, b, 
                            theta, epsilon2);
  }

  return err;
}
