#include "quadtree.h"

/* create a new quadtree */
quadtree_t *quadtree_new(double x1, double y1, 
                         double x2, double y2) {
  quadtree_t *quadtree = calloc(1, sizeof(quadtree_t));
  if (quadtree == NULL) return NULL;

  /* initialize lock for node */
  omp_init_lock(&(quadtree->node_lock));

  /* set quadtree bounding box */
  quadtree->x1 = x1;
  quadtree->y1 = y1;
  quadtree->x2 = x2;
  quadtree->y2 = y2;

  /* set quadtree center and side length */
  quadtree->l  = x2 - x1;
  quadtree->xm = 0.5 * (x1 + x2);
  quadtree->ym = 0.5 * (y1 + y2);

  /* set quadtree cumulative mass and center of mass */
  quadtree->work = 0;
  quadtree->m  = 0.0;
  quadtree->xc = 0.0;
  quadtree->yc = 0.0;

  /* initialize body at quadtree */
  quadtree->body = NULL;

  /* set quadtree subtrees */
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
int quadtree_insert(quadtree_t *quadtree, body_t *b) {
  double x = b->x, y = b->y;

  /* check if body is in the quadtree */
  if (!in_boundary(quadtree, x, y)) return 0;

  /* attempt to insert body at leaf */
  if (is_leaf(quadtree)) {
    
    /* immediately insert if quadtree has no tree */
    if (quadtree->body == NULL) {
      omp_set_lock(&(quadtree->node_lock));
      quadtree->body = b;
      omp_unset_lock(&(quadtree->node_lock));
      return 0;  
    } 

    /* there is a collision, we must allocate subtrees */
    else {
      double x1, y1, x2, y2, xm, ym, ox, oy;
      body_t *ob;

      /* get node values and remove node from tree */
      x1 = quadtree->x1;
      y1 = quadtree->y1;
      x2 = quadtree->x2;
      y2 = quadtree->y2;
      xm = quadtree->xm;
      ym = quadtree->ym;
      ob = quadtree->body;
      quadtree->body = NULL;
      ox = ob->x;
      oy = ob->y;

      /* initialize subtrees */
      omp_set_lock(&(quadtree->node_lock));
      quadtree->bot_left  = quadtree_new(x1, y1, xm, ym);
      quadtree->top_left  = quadtree_new(x1, ym, xm, y2);
      quadtree->bot_right = quadtree_new(xm, y1, x2, ym);
      quadtree->top_right = quadtree_new(xm, ym, x2, y2);
      omp_unset_lock(&(quadtree->node_lock));

      if (quadtree->bot_left  == NULL ||
          quadtree->top_left  == NULL ||
          quadtree->bot_right == NULL ||
          quadtree->top_right == NULL) return -1;

      /* re-insert removed node into subtree */
      quadtree_t *q1 = get_subtree(quadtree, ox, oy);
      if (q1 == NULL) return -1;
      if (quadtree_insert(q1, ob)) return -1;
    }
  }

  /* insert node into subtree */
  quadtree_t *q2 = get_subtree(quadtree, x, y);
  if (q2 == NULL) return -1;
  if (quadtree_insert(q2, b)) return -1;

  return 0;
}

/* post-order traversal of quadtree to compute 
 * node approximations and cumulative workloads */
int quadtree_traverse(quadtree_t *quadtree) {
  if (is_leaf(quadtree)) {
    body_t *b = quadtree->body;
    if (b != NULL) {
      quadtree->work = b->work;
      quadtree->m  = b->m;
      quadtree->xc = b->x;
      quadtree->yc = b->y;
    }
    return 0;
  }

  // int op;
  // #pragma omp parallel for schedule(static)
  // for (op = 0; op < 4; op++) {
  //   switch (op) {
  //     case 0:
  //       quadtree_traverse(quadtree->bot_left);
  //       break;
  //     case 1:
  //       quadtree_traverse(quadtree->top_left);
  //       break;
  //     case 2:
  //       quadtree_traverse(quadtree->bot_right);
  //       break;
  //     case 3:
  //     default:
  //       quadtree_traverse(quadtree->top_right);
  //       break;
  //   }
  // }

  quadtree_traverse(quadtree->bot_left);
  quadtree_traverse(quadtree->top_left);
  quadtree_traverse(quadtree->bot_right);
  quadtree_traverse(quadtree->top_right);

  double bot_left_m  = quadtree->bot_left->m,
         top_left_m  = quadtree->top_left->m,
         bot_right_m = quadtree->bot_right->m,
         top_right_m = quadtree->top_right->m, m, im;
  m = bot_left_m + top_left_m + bot_right_m + top_right_m;
  im = m > 1e-10 ? 1.0 / m : 1.0;

  quadtree->work = quadtree->bot_left->work  + 
                   quadtree->top_left->work +
                   quadtree->bot_right->work +
                   quadtree->top_right->work;
  quadtree->xc = (quadtree->bot_left->xc  * bot_left_m +
                  quadtree->top_left->xc  * top_left_m +
                  quadtree->bot_right->xc * bot_right_m +
                  quadtree->top_right->xc * top_right_m) * im;
  quadtree->yc = (quadtree->bot_left->yc  * bot_left_m +
                  quadtree->top_left->yc  * top_left_m +
                  quadtree->bot_right->yc * bot_right_m +
                  quadtree->top_right->yc * top_right_m) * im;
  quadtree->m = m;

  // #pragma omp parallel for schedule(static)
  // for (op = 0; op < 4; op++) {
  //   switch (op) {
  //     case 0:
  //       quadtree->work = quadtree->bot_left->work  + 
  //                        quadtree->top_left->work +
  //                        quadtree->bot_right->work +
  //                        quadtree->top_right->work;
  //       break;
  //     case 1:
  //       quadtree->xc = (quadtree->bot_left->xc  * bot_left_m +
  //                       quadtree->top_left->xc  * top_left_m +
  //                       quadtree->bot_right->xc * bot_right_m +
  //                       quadtree->top_right->xc * top_right_m) * im;
  //       break;
  //     case 2:
  //       quadtree->yc = (quadtree->bot_left->yc  * bot_left_m +
  //                       quadtree->top_left->yc  * top_left_m +
  //                       quadtree->bot_right->yc * bot_right_m +
  //                       quadtree->top_right->yc * top_right_m) * im;
  //       break;
  //     case 3:
  //     default:
  //       quadtree->m = m;
  //       break;
  //   }
  // }
  
  return 0;
}

/* post-order traversal of the tree to get the partition
 * of quadtree nodes for the given thread */
int quadtree_partition(quadtree_t *quadtree, partition_t *partition) {


  return 0;
}

/* aggregate forces for a body using the tree */
int quadtree_aggregate_forces(quadtree_t *quadtree, body_t *b,
                              double theta, double epsilon2) {
  double x = b->x, y = b->y, dx, dy, dmag,
         denom, denom2, mult;

  b->work++;

  /* check if the condition is satisfied */
  dx = quadtree->xc - x;
  dy = quadtree->yc - y;
  dmag = sqrt(dx*dx + dy*dy);

  if (quadtree->l/dmag < theta) {

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
  if (quadtree->body != NULL) {
    body_t *ob = quadtree->body;
    dx = ob->x - x;
    dy = ob->y - y;
    dmag = sqrt(dx*dx + dy*dy);

    denom = dmag*dmag + epsilon2;
    denom2 = denom*denom;
    denom = (denom * denom2) / denom2;
    mult = ob->m / denom;

    b->ax += mult * dx;
    b->ay += mult * dy;
  }

  /* if we are at a leaf, stop aggregating forces */
  if (is_leaf(quadtree)) return 0;

  /* recurse to aggregate subtrees */
  aggregate_forces(quadtree->bot_left, b, theta, epsilon2);
  aggregate_forces(quadtree->top_left, b, theta, epsilon2);
  aggregate_forces(quadtree->bot_right, b, theta, epsilon2);
  aggregate_forces(quadtree->top_right, b, theta, epsilon2);

  return 0;
}
