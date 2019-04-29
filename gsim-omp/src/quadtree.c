#include "quadtree.h"

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

/* free all of a quadtree's memory */
void quadtree_free(quadtree_t *quadtree) {
  if (is_leaf(quadtree)) { free(quadtree); return; }

  if (quadtree != NULL) {
    quadtree_free(quadtree->bot_left);
    quadtree_free(quadtree->top_left);
    quadtree_free(quadtree->bot_right);
    quadtree_free(quadtree->top_right);
    omp_destroy_lock(&(quadtree->node_lock));
    free(quadtree);
  }
}

/* insert a body into the quadtree */
void quadtree_insert(quadtree_t *quadtree, body_t *b) {
  double x = b->x, y = b->y;

  /* check if body is in the quadtree */
  if (!in_boundary(quadtree, x, y)) return;

  /* attempt to insert body at leaf */
  omp_set_lock(&(quadtree->node_lock));
  if (is_leaf(quadtree)) {
    
    /* immediately insert if quadtree has no tree */
    if (quadtree->body == NULL) {
      quadtree->body = b;
      omp_unset_lock(&(quadtree->node_lock));
      return;  
    } 
  
    /* there is a collision, we must allocate subtrees */
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
    quadtree->bot_left  = quadtree_new(x1, y1, xm, ym);
    quadtree->top_left  = quadtree_new(x1, ym, xm, y2);
    quadtree->bot_right = quadtree_new(xm, y1, x2, ym);
    quadtree->top_right = quadtree_new(xm, ym, x2, y2);
  
    omp_unset_lock(&(quadtree->node_lock));

    /* re-insert removed node into subtree */
    quadtree_t *q1 = get_subtree(quadtree, ox, oy);
    if (q1 == NULL) return;
    quadtree_insert(q1, ob);
  } else omp_unset_lock(&(quadtree->node_lock));

  /* insert node into subtree */
  quadtree_t *q2 = get_subtree(quadtree, x, y);
  if (q2 == NULL) return;
  quadtree_insert(q2, b);
}

/* helper function for quadtree_traverse with updating work */
static inline void combine_work_mass(quadtree_t *quadtree) {
  double bot_left_m  = quadtree->bot_left->m,
         top_left_m  = quadtree->top_left->m,
         bot_right_m = quadtree->bot_right->m,
         top_right_m = quadtree->top_right->m, m, im;
  m = bot_left_m + top_left_m + bot_right_m + top_right_m;
  im = m > 1e-10 ? 1.0 / m : 1.0;

  quadtree->work = quadtree->bot_left->work + 
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
}

/* post-order traversal of quadtree to compute 
 * node approximations and cumulative workloads */
void quadtree_traverse(quadtree_t *quadtree) {
  if (quadtree == NULL) return;
  if (is_leaf(quadtree)) {
    body_t *b = quadtree->body;
    if (b != NULL) {
      quadtree->work = b->work;
      quadtree->m  = b->m;
      quadtree->xc = b->x;
      quadtree->yc = b->y;
    }
    return;
  } 

  quadtree_traverse(quadtree->bot_left);
  quadtree_traverse(quadtree->top_left);
  quadtree_traverse(quadtree->bot_right);
  quadtree_traverse(quadtree->top_right);

  combine_work_mass(quadtree);
}

/* post-order traversal of the tree to get the partition
 * of quadtree nodes for the given thread */
void quadtree_partition(quadtree_t *quadtree,
                        partition_t *partition, int cur_work) {
  if (is_leaf(quadtree)) {
    body_t *b = quadtree->body;
    if (b != NULL && partition->min_work <= cur_work) {
      partition->pbodies[partition->num_pbodies++] = b;
    }
    return;
  }

  int min_work = partition->min_work,
      max_work = partition->max_work,
      next_work = 0;

  next_work = cur_work + quadtree->bot_left->work;
  if (min_work < next_work)
    quadtree_partition(quadtree->bot_left, partition, cur_work);
  if (max_work < next_work) return;
  cur_work = next_work;

  next_work = cur_work + quadtree->top_left->work;
  if (min_work < next_work)
    quadtree_partition(quadtree->top_left, partition, cur_work);
  if (max_work < next_work) return;
  cur_work = next_work;

  next_work = cur_work + quadtree->bot_right->work;
  if (min_work < next_work)
    quadtree_partition(quadtree->bot_right, partition, cur_work);
  if (max_work < next_work) return;
  cur_work = next_work;

  next_work = cur_work + quadtree->top_right->work;
  if (min_work < next_work)
    quadtree_partition(quadtree->top_right, partition, cur_work);
}

/* aggregate forces for a body using the tree */
void quadtree_aggregate_forces(quadtree_t *quadtree, body_t *b,
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

    #pragma omp atomic
    b->ax += mult * dx;
    #pragma omp atomic
    b->ay += mult * dy;

    return;
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

    #pragma omp atomic
    b->ax += mult * dx;
    #pragma omp atomic
    b->ay += mult * dy;
  }

  /* if we are at a leaf, stop aggregating forces */
  if (is_leaf(quadtree)) return;

  /* recurse to aggregate subtrees */
  quadtree_aggregate_forces(quadtree->bot_left, b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->top_left, b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->bot_right, b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->top_right, b, theta, epsilon2);
}
