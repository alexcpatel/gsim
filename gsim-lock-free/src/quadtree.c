#include "quadtree.h"

/* convert a body int to a pointer */
static inline body_t *int_to_body(uintptr_t body) {
  return (body_t *) (body & ~(uintptr_t) 0x1);
}

/* convert a body pointer to an int */
static inline uintptr_t body_to_int(body_t *b) {
  return (uintptr_t) b | (uintptr_t) 0x1;
}

/* get the next subtree for the position */
static inline quadtree_t *get_subtree(quadtree_t *quadtree,
                                      double x, double y) {
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
  if (!quadtree) return NULL;

  /* initialize body at quadtree */
  quadtree->body = 0x1; // is leaf node, NULL body

  /* set quadtree bounding box */
  quadtree->x1 = x1;
  quadtree->y1 = y1;
  quadtree->x2 = x2;
  quadtree->y2 = y2;

  /* set quadtree center and side length */
  quadtree->l  = x2 - x1;
  quadtree->xm = 0.5 * (x1 + x2);
  quadtree->ym = 0.5 * (y1 + y2);

  /* set quadtree subtrees */
  quadtree->top_left  = NULL;
  quadtree->top_right = NULL;
  quadtree->bot_left  = NULL;
  quadtree->bot_right = NULL;

  /* set accumulated quadtree info */
  quadtree->m  = 0.0;
  quadtree->xc = 0.0;
  quadtree->yc = 0.0;
  quadtree->work = 0;

  return quadtree;
}

/* free all of a quadtree's memory */
void quadtree_free(quadtree_t *quadtree) {
  if (quadtree) {
    if (!(quadtree->body & (uintptr_t) 0x1)) {
      quadtree_free(quadtree->bot_left);
      quadtree_free(quadtree->top_left);
      quadtree_free(quadtree->bot_right);
      quadtree_free(quadtree->top_right);
    }
    free(quadtree);
  }
}

/* insert a body into the quadtree */
void quadtree_insert(quadtree_t *quadtree, body_t *b) {
  double x = b->x, y = b->y;

  /* check if body is in the quadtree */
  if (quadtree == NULL ||
      x < quadtree->x1 ||
      x > quadtree->x2 ||
      y < quadtree->y1 ||
      y > quadtree->y2) return;
  
  /* if this leaf is empty, insert immediately */
  while (quadtree->body == (uintptr_t) 0x1) {
    uintptr_t nib = body_to_int(b);
    if (__sync_bool_compare_and_swap(&(quadtree->body), 
        (uintptr_t) 0x1, nib)) return;
  }

  /* if there is a collision, allocate subtrees */
  uintptr_t oib;
  while ((oib = quadtree->body)   &&
         (oib & ~(uintptr_t) 0x1) &&
          oib &  (uintptr_t) 0x1) {
    double x1, y1, x2, y2, xm, ym;

    /* get node values and remove node from tree */
    x1 = quadtree->x1;
    y1 = quadtree->y1;
    x2 = quadtree->x2;
    y2 = quadtree->y2;
    xm = quadtree->xm;
    ym = quadtree->ym;

    /* initialize subtrees */
    if (!quadtree->bot_left) {
      quadtree_t *bot_left = quadtree_new(x1, y1, xm, ym);
      while (1) {
        if (__sync_bool_compare_and_swap(&(quadtree->bot_left), 
            NULL, bot_left)) break;
        if (quadtree->bot_left) {
          free(bot_left);
          break;
        }
      }
    }

    if (!quadtree->top_left) {
      quadtree_t *top_left = quadtree_new(x1, ym, xm, y2);
      while (1) {
        if (__sync_bool_compare_and_swap(&(quadtree->top_left), 
            NULL, top_left)) break;
        if (quadtree->top_left) {
          free(top_left);
          break;
        }
      }
    }

    if (!quadtree->bot_right) {
      quadtree_t *bot_right = quadtree_new(xm, y1, x2, ym);
      while (1) {
        if (__sync_bool_compare_and_swap(&(quadtree->bot_right), 
            NULL, bot_right)) break;
        if (quadtree->bot_right) {
          free(bot_right);
          break;
        }
      }
    }

    if (!quadtree->top_right) {
      quadtree_t *top_right = quadtree_new(xm, ym, x2, y2);
      while (1) {
        if (__sync_bool_compare_and_swap(&(quadtree->top_right), 
            NULL, top_right)) break;
        if (quadtree->top_right) {
          free(top_right);
          break;
        }
      }
    }

    /* get body to re-insert */
    body_t *ob = int_to_body(oib);
    if (__sync_bool_compare_and_swap(&(quadtree->body), 
        oib, (uintptr_t) 0x0)) {
      quadtree_insert(get_subtree(quadtree, ob->x, ob->y), ob);
    }
  }

  /* insert node into subtree */
  quadtree_insert(get_subtree(quadtree, x, y), b);
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
  if (!quadtree) return;
  if (quadtree->body & (uintptr_t) 0x1) {
    body_t *b = int_to_body(quadtree->body);
    if (b) {
      quadtree->m  = b->m;
      quadtree->xc = b->x;
      quadtree->yc = b->y;
      quadtree->work = b->work;
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
  if (quadtree->body & (uintptr_t) 0x1) {
    body_t *b = int_to_body(quadtree->body);
    if (b && partition->min_work <= cur_work) {
      partition->pbodies[partition->num_pbodies++] = b;
    }
    return;
  }

  uint32_t min_work = partition->min_work,
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
  double x = b->x, y = b->y, dx, dy, dmag;
  b->work++;

  /* check if the condition is satisfied */
  dx = quadtree->xc - x;
  dy = quadtree->yc - y;
  dmag = sqrt(dx*dx + dy*dy);

  /* add subtree contribution */
  if (quadtree->l/dmag < theta) {
    double mult = quadtree->m / (dmag*dmag + epsilon2);

    b->ax += mult * dx;
    b->ay += mult * dy;
    return;
  }

  /* add node contribution */
  body_t *ob = int_to_body(quadtree->body);
  if (ob) {
    dx = ob->x - x;
    dy = ob->y - y;
    dmag = sqrt(dx*dx + dy*dy);

    double mult = ob->m / (dmag*dmag + epsilon2);

    b->ax += mult * dx;
    b->ay += mult * dy;
  }

  /* if we are at a leaf, stop aggregating forces */
  if (quadtree->body & (uintptr_t) 0x1) return;

  /* recurse to aggregate subtrees */
  quadtree_aggregate_forces(quadtree->bot_left,  b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->top_left,  b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->bot_right, b, theta, epsilon2);
  quadtree_aggregate_forces(quadtree->top_right, b, theta, epsilon2);
}
