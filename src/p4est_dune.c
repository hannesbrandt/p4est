/*
  This file is part of p4est.
  p4est is a C library to manage a collection (a forest) of multiple
  connected adaptive quadtrees or octrees in parallel.

  Copyright (C) 2010 The University of Texas System
  Additional copyright (C) 2011 individual authors
  Written by Carsten Burstedde, Lucas C. Wilcox, and Tobin Isaac

  p4est is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  p4est is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with p4est; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef P4_TO_P8
#include <p4est_algorithms.h>
#include <p4est_bits.h>
#include <p4est_dune.h>
#include <p4est_lnodes.h>
#include <p4est_search.h>
#else
#include <p8est_algorithms.h>
#include <p8est_bits.h>
#include <p8est_dune.h>
#include <p8est_lnodes.h>
#include <p8est_search.h>
#endif

/* map element corner to element node index */
static const int    corner_enode[P4EST_CHILDREN] = {
  0, 4, 20, 24
#ifdef P4_TO_P8
    , 100, 104, 120, 124
#endif
};

#ifdef P4_TO_P8

/* map element edge to element node index */
static const int    edge_enode[P8EST_EDGES] = {
  2, 22, 102, 122, 10, 14, 110, 114, 50, 54, 70, 74,
};

#endif

/* map element face to element node index */
static const int    face_enode[P4EST_FACES] = {
#ifndef P4_TO_P8
  10, 14, 2, 22
#else
  60, 64, 52, 72, 12, 112
#endif
};

/* element node next to a corner c one step along a face of given normal i */
static int
corner_face_enode (int c, int i)
{
  int                 j;
  int                 m;
  int                 eno;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= i && i < P4EST_DIM);

  /* move by one step in each tangential direction */
  eno = corner_enode[c];
  m = 1;
  for (j = 0; j < P4EST_DIM; ++j) {
    if (j != i) {
      eno += (1 - 2 * (c & 1)) * m;
    }
    m *= 5;
    c >>= 1;
  }
  P4EST_ASSERT (c == 0);

  /* return final element node number */
  P4EST_ASSERT (0 <= eno && eno < m);
  return eno;
}

#ifdef P4_TO_P8

/* element node next to a corner c along j and two steps along k */
static int
corner_twostep_enode (int c, int j, int k)
{
  int                 i;
  int                 m;
  int                 eno;

  P4EST_ASSERT (0 <= c && c < P4EST_CHILDREN);
  P4EST_ASSERT (0 <= j && j < P4EST_DIM);
  P4EST_ASSERT ((0 <= k && k < P4EST_DIM) || k == -1);
  P4EST_ASSERT (j != k);

  /* move by one or two steps depending on direction */
  eno = corner_enode[c];
  m = 1;
  for (i = 0; i < P4EST_DIM; ++i) {
    if (i == j) {
      eno += (1 - 2 * (c & 1)) * m;
    }
    else if (i == k) {
      eno += (1 - 2 * (c & 1)) * 2 * m;
    }
    m *= 5;
    c >>= 1;
  }
  P4EST_ASSERT (c == 0);

  /* return final element node number */
  P4EST_ASSERT (0 <= eno && eno < m);
  return eno;
}

/* element node next to a corner c one step along a given edge j */
static int
corner_edge_enode (int c, int j)
{
  return corner_twostep_enode (c, j, -1);
}

#endif /* P4_TO_P8 */

/* set default parameters */
void
p4est_dune_numbers_params_init (p4est_dune_numbers_params_t * params)
{
  P4EST_ASSERT (params != NULL);

  /* set default parameters to construct a dune numbers structure */
  params->ctype = P4EST_CONNECT_FULL;
}

#if 0                           /* presently unused */

static              p4est_gloidx_t
lni_to_gni (p4est_lnodes_t * ln, p4est_locidx_t lni)
{
  p4est_locidx_t      oco;

  P4EST_ASSERT (0 <= lni);

  if (lni < (oco = ln->owned_count)) {
    return ln->global_offset + lni;
  }
  else {
    P4EST_ASSERT (lni < ln->num_local_nodes);
    return ln->nonlocal_nodes[lni - oco];
  }
}

#endif

/* Populate number arrays for real.  This function does all the work. */
static void
generate_numbers (p4est_dune_numbers_t * dn, p4est_lnodes_t * ln)
{
  p4est_lnodes_code_t fc;
  p4est_locidx_t      lne, el;
  p4est_locidx_t     *enos;
  int                 i;
  int                 cid, work;
  int                 c, f;
  int                 do_c, do_f;
  p4est_locidx_t     *geco, *gefa;
#ifdef P4_TO_P8
  int                 j, k;
  int                 e;
  int                 do_e;
  p4est_locidx_t     *geed;
#endif

  /* basic verification of input parameters */
  P4EST_ASSERT (dn != NULL);
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->degree == 4);

  /* examine which codimensions are handled */
  do_c = dn->element_corners != NULL;
  geco = NULL;
#ifdef P4_TO_P8
  do_e = dn->element_edges != NULL;
  geed = NULL;
#endif
  do_f = dn->element_faces != NULL;
  gefa = NULL;

  /* maybe there is nothing to do */
  if (!do_c &&
#ifdef P4_TO_P8
      !do_e &&
#endif
      !do_f) {
    /* there is absolutely nothing to do */
    return;
  }

  /* prepare loop over elements to assign globally unique node indices */
  lne = ln->num_local_elements;
  enos = ln->element_nodes;
  for (el = 0; el < lne; ++el) {

    /* assign standard values for an element that is not hanging anywhere */
    if (do_c) {
      geco = (p4est_locidx_t *)
        sc_array_index (dn->element_corners, el * P4EST_CHILDREN);
      for (c = 0; c < P4EST_CHILDREN; ++c) {
        /* retrieve global number of corner node */
        geco[c] = enos[corner_enode[c]];
      }
    }
#ifdef P4_TO_P8
    if (do_e) {
      geed = (p4est_locidx_t *)
        sc_array_index (dn->element_edges, el * P8EST_EDGES);
      for (e = 0; e < P8EST_EDGES; ++e) {
        /* retrieve global number of edge node */
        geed[e] = enos[edge_enode[e]];
      }
    }
#endif
    if (do_f) {
      gefa = (p4est_locidx_t *)
        sc_array_index (dn->element_faces, el * P4EST_FACES);
      for (f = 0; f < P4EST_FACES; ++f) {
        /* retrieve global number of face node */
        gefa[f] = enos[face_enode[f]];
      }
    }

    /* modify node numbers for hanging faces and edges */
    if (0 != (fc = ln->face_code[el])) {
      /* element has child number cid and some hanging faces or edges */
      cid = (fc & (P4EST_CHILDREN - 1));
      work = (fc >> P4EST_DIM);

      /* iterate over the faces adjacent to corner cid */
      for (i = 0; i < P4EST_DIM; ++i) {
        if (work & 1) {
          /* face normal to direction i is hanging */
          f = p4est_corner_faces[cid][i];
          c = cid ^ (P4EST_CHILDREN - 1) ^ (1 << i);
          if (do_c) {
            /* assign to mid-face corner */
            geco[c] = enos[face_enode[f]];
          }
#ifdef P4_TO_P8
          if (do_e) {
            /* assign to the two mid-face edges */
            j = (i + 1) % 3;
            k = (j + 1) % 3;
            geed[p8est_corner_edges[c][j]] =
              enos[corner_twostep_enode (cid, j, k)];
            geed[p8est_corner_edges[c][k]] =
              enos[corner_twostep_enode (cid, k, j)];
          }
#endif
          if (do_f) {
            /* assign to small mid-face */
            gefa[f] = enos[corner_face_enode (cid, i)];
          }
        }
        work >>= 1;
      }
#ifdef P4_TO_P8
      /* iterate over the edges adjacent to corner cid */
      for (j = 0; j < P4EST_DIM; ++j) {
        if (work & 1) {
          /* edge in direction j is hanging */
          e = p8est_corner_edges[cid][j];
          if (do_c) {
            /* assign to mid-edge corner */
            c = cid ^ (1 << j);
            geco[c] = enos[edge_enode[e]];
          }
          if (do_e) {
            /* assign to small mid-edge */
            geed[e] = enos[corner_edge_enode (cid, j)];
          }
        }
        work >>= 1;
      }
#endif
      P4EST_ASSERT (work == 0);
    }

    /* proceed to the next element */
    enos += ln->vnodes;
  }
}

static void
consecutive_numbers (p4est_dune_numbers_t * dn,
                     p4est_locidx_t num_local_nodes)
{
  /* an array of bits, one per local node number, for codimensions */
  uint8_t            *local_node_bytes;
  uint8_t             mask_byte, node_byte;
  int                 j;
  sc_array_t         *dim_numbers[P4EST_DIM], *numbers;
  p4est_locidx_t     *dim_results[P4EST_DIM], *newnums;
  p4est_locidx_t      li, lcount, lnum, ncount;

  /* input checks */
  P4EST_ASSERT (dn != NULL);

  /* prepare temporary data */
  local_node_bytes = P4EST_ALLOC_ZERO (uint8_t, num_local_nodes);
  dim_numbers[0] = dn->element_corners;
  dim_results[0] = &dn->num_corner_numbers;
#ifdef P4_TO_P8
  dim_numbers[1] = dn->element_edges;
  dim_results[1] = &dn->num_edge_numbers;
#endif
  dim_numbers[P4EST_DIM - 1] = dn->element_faces;
  dim_results[P4EST_DIM - 1] = &dn->num_face_numbers;

  /* go through node arrays to count unique entries */
  for (j = 0; j < P4EST_DIM; ++j) {
    P4EST_ASSERT (*dim_results[j] == 0);
    if ((numbers = dim_numbers[j]) == NULL) {
      continue;
    }
    P4EST_ASSERT (numbers->elem_size == sizeof (p4est_locidx_t));
    mask_byte = 1 << j;

    /* mark each local node that occurs in any of the arrays */
    lcount = (p4est_locidx_t) numbers->elem_count;
    for (li = 0; li < lcount; ++li) {
      lnum = *(p4est_locidx_t *) sc_array_index (numbers, li);
      P4EST_ASSERT (0 <= lnum && lnum < num_local_nodes);
      node_byte = local_node_bytes[lnum];
      if (!(node_byte & mask_byte)) {
        local_node_bytes[lnum] = node_byte | mask_byte;
      }
    }
  }

  /* renumber nodes to contiguous subrange */
  newnums = P4EST_ALLOC (p4est_locidx_t, num_local_nodes);
  for (j = 0; j < P4EST_DIM; ++j) {
    P4EST_ASSERT (*dim_results[j] == 0);
    if ((numbers = dim_numbers[j]) == NULL) {
      continue;
    }
    mask_byte = 1 << j;
    memset (newnums, -1, num_local_nodes * sizeof (*newnums));
    ncount = 0;

    /* figure out which local nodes are relevant for dimension j */
    lcount = num_local_nodes;
    for (li = 0; li < lcount; ++li) {
      if (local_node_bytes[li] & mask_byte) {
        newnums[li] = ncount++;
      }
    }

    /* reassign entries of numbers array */
    lcount = (p4est_locidx_t) numbers->elem_count;
    for (li = 0; li < lcount; ++li) {
      lnum = *(p4est_locidx_t *) sc_array_index (numbers, li);
      P4EST_ASSERT (0 <= lnum && lnum < num_local_nodes);
      P4EST_ASSERT (local_node_bytes[lnum] & mask_byte);
      lnum = newnums[lnum];
      P4EST_ASSERT (0 <= lnum && lnum < ncount);
      *(p4est_locidx_t *) sc_array_index (numbers, li) = lnum;
    }
    *dim_results[j] = ncount;
  }
  P4EST_FREE (newnums);
  P4EST_FREE (local_node_bytes);
}

/* create lookup tables for unique corners and faces */
p4est_dune_numbers_t *
p4est_dune_numbers_new (p4est_t * p4est, p4est_ghost_t * ghost,
                        const p4est_dune_numbers_params_t * params)
{
  p4est_dune_numbers_t *dn;
  p4est_dune_numbers_params_t *pa;
  p4est_lnodes_t     *ln;
  p4est_locidx_t      lne, lnums;

  /* formally verify arguments */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est_is_valid (p4est));
  P4EST_ASSERT (p4est_is_balanced (p4est, P4EST_CONNECT_ALMOST));
  P4EST_ASSERT (ghost != NULL);

  /* allocate numbers structure to populate */
  dn = P4EST_ALLOC_ZERO (p4est_dune_numbers_t, 1);
  pa = &dn->params;

  /* establish input parameters and remember them */
  if (params == NULL) {
    p4est_dune_numbers_params_init (pa);
  }
  else {
    *pa = *params;
    params = NULL;
  }

  /* generate temporary local ghost (if needed) and node numbering */
  ln = p4est_lnodes_new (p4est, ghost, 4);

  /* construct corner and face numbers and return them */
  lne = p4est->local_num_quadrants;
  if (pa->ctype >= P4EST_CONNECT_CORNER) {
    dn->element_corners =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P4EST_CHILDREN);
  }
#ifdef P4_TO_P8
  if (pa->ctype >= P8EST_CONNECT_EDGE) {
    dn->element_edges =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P8EST_EDGES);
  }
#endif
  if (pa->ctype >= P4EST_CONNECT_FACE) {
    dn->element_faces =
      sc_array_new_count (sizeof (p4est_locidx_t), lne * P4EST_FACES);
  }
  generate_numbers (dn, ln);

  /* free temporary data */
  lnums = ln->num_local_nodes;
  p4est_lnodes_destroy (ln);

  /* transform numbering into contiguous ranges per codimension */
  consecutive_numbers (dn, lnums);

  /* that's it! */
  return dn;
}

void
p4est_dune_numbers_destroy (p4est_dune_numbers_t * dn)
{
  P4EST_ASSERT (dn != NULL);

  if (dn->element_corners != NULL) {
    sc_array_destroy (dn->element_corners);
  }
#ifdef P4_TO_P8
  if (dn->element_edges != NULL) {
    sc_array_destroy (dn->element_edges);
  }
#endif
  if (dn->element_faces != NULL) {
    sc_array_destroy (dn->element_faces);
  }
  P4EST_FREE (dn);
}

typedef struct p4est_dune_iter
{
  /* wrapped iteration context */
  p4est_t            *p4est;
  p4est_iter_volume_t iter_volume;
  p4est_iter_face_t   iter_face;
  void               *user_data;

  /* this structure will be reused many times */
  p4est_iter_face_info_t sfinfo;
  p4est_iter_face_info_t *finfo;
}
p4est_dune_iter_t;

static void
dune_iter_volume (p4est_iter_volume_info_t * info, void *user_data)
{
  p4est_dune_iter_t  *dune_iter = (p4est_dune_iter_t *) user_data;

  /* consistency checks */
  P4EST_ASSERT (info != NULL);
  P4EST_ASSERT (dune_iter != NULL);
  P4EST_ASSERT (dune_iter->p4est == info->p4est);

  /* unwrap volume iteration callback */
  P4EST_ASSERT (dune_iter->iter_volume != NULL);
  dune_iter->iter_volume (info, dune_iter->user_data);
}

static void
dune_iter_face (p4est_iter_face_info_t * info, void *user_data)
{
  p4est_dune_iter_t  *dune_iter = (p4est_dune_iter_t *) user_data;
  p4est_iter_face_info_t *finfo;
  p4est_iter_face_side_t *s[2], *sh;
  p4est_iter_face_side_t *f[2], *fh;
  int                 i;
  int                 fside, hside;

  /* consistency checks */
  P4EST_ASSERT (info != NULL);
  P4EST_ASSERT (dune_iter != NULL);
  P4EST_ASSERT (dune_iter->p4est == info->p4est);

  /* unwrap face callback */
  P4EST_ASSERT (dune_iter->iter_face != NULL);

  /* examine face boundary situation */
  P4EST_ASSERT (info->sides.elem_count > 0);
  if (info->sides.elem_count == 1) {
    /* this face is on a physical boundary */
    dune_iter->iter_face (info, dune_iter->user_data);
    return;
  }

  /* examine hanging face situation */
  P4EST_ASSERT (info->sides.elem_count == 2);
  for (i = 0; i < 2; ++i) {
    s[i] = (p4est_iter_face_side_t *) sc_array_index_int (&info->sides, i);
    P4EST_ASSERT (0 <= s[i]->face && s[i]->face < P4EST_FACES);
  }
  if (!s[0]->is_hanging && !s[1]->is_hanging) {
    dune_iter->iter_face (info, dune_iter->user_data);
    return;
  }

  /* we have precisely one hanging face.  Reconstruct its info */
  P4EST_ASSERT (s[0]->is_hanging != s[1]->is_hanging);
  finfo = dune_iter->finfo;
  P4EST_ASSERT (finfo != NULL);
  P4EST_ASSERT (finfo->sides.elem_size == sizeof (p4est_iter_face_side_t));
  P4EST_ASSERT (finfo->sides.elem_count == 2);
  if (s[0]->is_hanging) {
    hside = 0;
    fside = 1;
  }
  else {
    hside = 1;
    fside = 0;
  }
  P4EST_ASSERT (s[hside]->is_hanging);
  P4EST_ASSERT (!s[fside]->is_hanging);
  for (i = 0; i < 2; ++i) {
    f[i] = (p4est_iter_face_side_t *) sc_array_index_int (&finfo->sides, i);
  }

  /* copy full side information as is and prepare setting the other side */
  *f[fside] = *s[fside];
  sh = s[hside];
  fh = f[hside];
  fh->treeid = sh->treeid;
  fh->face = sh->face;
  fh->is_hanging = sh->is_hanging;
  P4EST_ASSERT (fh->is_hanging);

  /* iterate over the hanging face neighbors */
  for (i = 0; i < P4EST_HALF; ++i) {
    if (s[fside]->is.full.is_ghost && sh->is.hanging.is_ghost[i]) {
      /* if both sides are remote, we do nothing for this pair */
      continue;
    }

    /* we copy the hanging face into the first position */
    fh->is.hanging.is_ghost[0] = sh->is.hanging.is_ghost[i];
    fh->is.hanging.quad[0] = sh->is.hanging.quad[i];
    fh->is.hanging.quadid[0] = sh->is.hanging.quadid[i];

    /* abuse the second position to store this face's sequence number */
    fh->is.hanging.quadid[1] = i;

    /* execute callback for this face pairing */
    dune_iter->iter_face (finfo, dune_iter->user_data);
  }
}

static void
p4est_dune_iterate_wrap (p4est_t * p4est, p4est_ghost_t * ghost_layer,
                         void *user_data, p4est_iter_volume_t iter_volume,
                         p4est_iter_face_t iter_face)
{
  P4EST_ASSERT (p4est != NULL);

  if (iter_face == NULL) {
    /* no need to do anything special: pass through */
    p4est_iterate (p4est, ghost_layer, user_data, iter_volume, NULL,
#ifdef P4_TO_P8
                   NULL,
#endif
                   NULL);
  }
  else {
    /* we need to translate the face iteration calls */
    p4est_iter_face_info_t *finfo;
    p4est_dune_iter_t   sdune_iter, *dune_iter = &sdune_iter;
    memset (dune_iter, 0, sizeof (*dune_iter));

    /* remember wrapped information */
    dune_iter->p4est = p4est;
    dune_iter->iter_volume = iter_volume;
    dune_iter->iter_face = iter_face;
    dune_iter->user_data = user_data;

    /* prepare internal context */
    finfo = dune_iter->finfo = &dune_iter->sfinfo;
    finfo->p4est = p4est;
    finfo->ghost_layer = ghost_layer;
    sc_array_init_count (&finfo->sides, sizeof (p4est_iter_face_side_t), 2);
    sc_array_memset (&finfo->sides, 0);

    /* wrapped iterator call */
    p4est_iterate (p4est, ghost_layer, dune_iter,
                   iter_volume != NULL ? dune_iter_volume : NULL,
                   /* we know this is not NULL here, regardless */
                   iter_face != NULL ? dune_iter_face : NULL,
#ifdef P4_TO_P8
                   NULL,
#endif
                   NULL);

    /* free work memory */
    sc_array_reset (&finfo->sides);
  }
}

/**********************************************************************\
*                         Auxiliary hash cache                         *
\**********************************************************************/

typedef struct sc_dlink
{
  void               *data;
  struct sc_dlink    *prev;
  struct sc_dlink    *next;
}
sc_dlink_t;

typedef void        (*sc_drop_function_t) (void *v, void *u);

typedef struct sc_hash_mru
{
  /* functions provided by the user */
  sc_hash_function_t  hash_fn;
  sc_equal_function_t equal_fn;
  sc_drop_function_t  drop_fn;
  void               *user;

  /* internal container objects */
  sc_hash_t          *hash;
  sc_mempool_t       *pool;
  sc_dlink_t         *first, *last;

  /* counters and statistics */
  size_t              maxcount;
  size_t              count;
}
sc_hash_mru_t;

sc_hash_mru_t      *sc_hash_mru_new (sc_hash_function_t hash_fn,
                                     sc_equal_function_t equal_fn,
                                     sc_drop_function_t drop_fn,
                                     void *user, size_t maxcount);
void                sc_hash_mru_destroy (sc_hash_mru_t *mru);

int                 sc_hash_mru_insert_unique (sc_hash_mru_t *mru,
                                               void *v, void ***found);
int                 sc_hash_mru_remove (sc_hash_mru_t *mru,
                                        void *v, void **found);

#ifndef P4_TO_P8

static unsigned int
sc_hash_mru_hash (const void *v, const void *u)
{
  const sc_hash_mru_t *mru = (const sc_hash_mru_t *) u;
  const sc_dlink_t   *lynk = (const sc_dlink_t *) v;

  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->hash_fn != NULL);

  return mru->hash_fn (lynk->data, mru->user);
}

static int
sc_hash_mru_is_equal (const void *v1, const void *v2, const void *u)
{
  const sc_hash_mru_t *mru = (const sc_hash_mru_t *) u;
  const sc_dlink_t   *lynk1 = (const sc_dlink_t *) v1;
  const sc_dlink_t   *lynk2 = (const sc_dlink_t *) v2;

  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->equal_fn != NULL);

  return mru->equal_fn (lynk1->data, lynk2->data, mru->user);
}

static void
sc_hash_mru_consolidate (sc_hash_mru_t *mru)
{
  sc_dlink_t         *drop;

  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* drop superfluous objects */
  while (mru->count > mru->maxcount) {
    drop = mru->first;
    SC_ASSERT (drop != NULL);
    SC_ASSERT (drop->prev == NULL);

    /* call the user's drop handler */
    if (mru->drop_fn != NULL) {
      mru->drop_fn (drop->data, mru->user);
    }

    /* drop oldest list entry */
    if ((mru->first = drop->next) == NULL) {
      SC_ASSERT (mru->count == 1);
      mru->last = NULL;
    }
    else {
      SC_ASSERT (drop->next->prev == drop);
      mru->first->prev = NULL;
    }

    /* update memory and count */
    sc_mempool_free (mru->pool, drop);
    --mru->count;
  }
}

sc_hash_mru_t      *
sc_hash_mru_new (sc_hash_function_t hash_fn, sc_equal_function_t equal_fn,
                 sc_drop_function_t drop_fn, void *user, size_t maxcount)
{
  sc_hash_mru_t      *mru;

  SC_ASSERT (hash_fn != NULL);
  SC_ASSERT (equal_fn != NULL);

  mru = SC_ALLOC_ZERO (sc_hash_mru_t, 1);
  mru->hash_fn = hash_fn;
  mru->equal_fn = equal_fn;
  mru->drop_fn = drop_fn;
  mru->user = user;

  mru->hash = sc_hash_new (sc_hash_mru_hash, sc_hash_mru_is_equal, mru, NULL);
  mru->pool = sc_mempool_new (sizeof (sc_dlink_t));

  mru->maxcount = maxcount;

  return mru;
}

void
sc_hash_mru_destroy (sc_hash_mru_t *mru)
{
  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* call drop handler on remaining items */
  if (mru->drop_fn != NULL) {
    sc_dlink_t         *head = mru->first;

    /* walk through the list from oldest to newest */
    while (head != NULL) {
      mru->drop_fn (head->data, mru->user);
      head = head->next;

      /* returning to mempool would be redundant here */
    }
  }

  /* free all stored list elements */
  sc_hash_destroy (mru->hash);

  /* free the hash structure itself */
  sc_mempool_destroy (mru->pool);

  /* free this object */
  SC_FREE (mru);
}

int
sc_hash_mru_insert_unique (sc_hash_mru_t *mru, void *v, void ***found)
{
  int                 inserted;
  void              **lfound;
  sc_dlink_t          key, *lkey = &key;
  sc_dlink_t         *add;

  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* construct hash key */
  lkey->data = v;
  inserted = sc_hash_insert_unique (mru->hash, lkey, &lfound);
  if (inserted) {

    /* this object is newly added */
    add = (sc_dlink_t *) sc_mempool_alloc (mru->pool);
    add->data = v;
    add->next = NULL;
    if (mru->last == NULL) {

      /* the list was empty before */
      SC_ASSERT (mru->first == NULL && mru->count == 0);
      (mru->first = add)->prev = NULL;
    }
    else {

      /* append to the list */
      SC_ASSERT (mru->last->next == NULL && mru->count > 0);
      (mru->last->next = add)->prev = mru->last;
    }
    *(sc_dlink_t **) lfound = mru->last = add;
    ++mru->count;
  }
  else {

    /* this object exists already */
    add = *(sc_dlink_t **) lfound;
    if (add != mru->last) {

      /* remove it from its place */
      SC_ASSERT (add->next != NULL);
      if ((add->next->prev = add->prev) == NULL) {

        /* we are removing the first element */
        SC_ASSERT (add == mru->first);
        mru->first = add->next;
      }
      else {

        /* we keep the first element */
        add->prev->next = add->next;
      }

      /* and append it to the end */
      (add->prev = mru->last)->next = add;
      (mru->last = add)->next = NULL;
    }
  }

  /* return data location if so desired */
  if (found != NULL) {
    *found = &add->data;
  }

  /* indicate pre-existing object and return */
  sc_hash_mru_consolidate (mru);
  return inserted;
}

/* not calling drop on the item found! */
int
sc_hash_mru_remove (sc_hash_mru_t *mru, void *v, void **found)
{
  int                 removed;
  void               *lfound;
  sc_dlink_t          key, *lkey = &key;
  sc_dlink_t         *drop;

  /* verify preconditions */
  SC_ASSERT (mru != NULL);
  SC_ASSERT (mru->pool->elem_count == mru->count);
  SC_ASSERT (mru->hash->elem_count == mru->count);

  /* construct hash key */
  lkey->data = v;
  removed = sc_hash_remove (mru->hash, lkey, &lfound);
  if (removed) {
    SC_ASSERT (mru->count > 0);

    /* return data location if so desired */
    drop = (sc_dlink_t *) lfound;
    if (found != NULL) {
      *found = drop->data;
    }

    /* unlink found object from list */
    if (drop->prev == NULL) {
      SC_ASSERT (drop == mru->first);
      mru->first = drop->next;
    }
    else {
      drop->prev->next = drop->next;
    }
    if (drop->next == NULL) {
      SC_ASSERT (drop == mru->last);
      mru->last = drop->prev;
    }
    else {
      drop->next->prev = drop->prev;
    }

    /* update memory and count */
    sc_mempool_free (mru->pool, drop);
    --mru->count;
  }

  /* indicate pre-existing object and return */
  sc_hash_mru_consolidate (mru);
  return removed;
}

#endif

/**********************************************************************\
*                       Non-balanced face iterator                     *
\**********************************************************************/

/** The complete quadrant context for the recursion. */
typedef struct p4est_quad_nonb
{
  /* Valid quadrant with valid p.which_tree member. */
  p4est_quadrant_t    skey;

  /* Number relative to tree of first local descendant of skey. */
  p4est_locidx_t      quadid;

  /* View on local quadrants that are descendants of skey. */
  sc_array_t          squads;

  /** If skey is equal to a leaf, this is all zeroes. */
  size_t              split[P4EST_CHILDREN + 1];

  /* Number of first ghost that is a descendant of skey. */
  p4est_locidx_t      ghostid;

  /* View on ghosts that are descendants of skey. */
  sc_array_t          ghosts;

  /* If skey is equal to a ghost, this is all zeroes. */
  size_t              gsplit[P4EST_CHILDREN + 1];

  /* Number of visible strict descendants below skey.
   * If this is 0, there is neither a local nor ghost contained.
   * The special value -1 designates a full size local quadrant.
   * The special value -2 designates a full size ghost quadrant.
   */
  p4est_locidx_t      nvdesc;
}
p4est_quad_nonb_t;

static unsigned
p4est_quad_nonb_hash (const void *v, const void *u)
{
  const p4est_quad_nonb_t *n = (const p4est_quad_nonb_t *) v;

  P4EST_ASSERT (n != NULL);

  return p4est_quadrant_hash_piggy (&n->skey);
}

static int
p4est_quad_nonb_is_equal (const void *v1, const void *v2, const void *u)
{
  const p4est_quad_nonb_t *n1 = (const p4est_quad_nonb_t *) v1;
  const p4est_quad_nonb_t *n2 = (const p4est_quad_nonb_t *) v2;

  P4EST_ASSERT (n1 != NULL && n2 != NULL);

  return p4est_quadrant_is_equal_piggy (&n1->skey, &n2->skey);
}

/** Global context object for non-balanced volume & face iterator. */
typedef struct p4est_dune_nonb
{
  /* general context */
  p4est_t            *p4est;
  p4est_ghost_t      *ghost_layer;
  void               *user_data;
  p4est_iter_volume_t iter_volume;
  p4est_iter_face_t   iter_face;

  /* containers and data */
  sc_hash_t          *qhash;
  sc_mempool_t       *qpool;

  /* newly developed hash cache */
  sc_hash_mru_t      *mru[P4EST_MAXLEVEL];

  /* state of recursion */
  p4est_tree_t       *tree;
  p4est_iter_volume_info_t vinfo;
  p4est_iter_face_info_t finfo;
}
p4est_dune_nonb_t;

static p4est_quad_nonb_t *
p4est_dune_nquad_alloc (p4est_dune_nonb_t *nonb)
{
  P4EST_ASSERT (nonb != NULL && nonb->qpool != NULL);
  return (p4est_quad_nonb_t *) sc_mempool_alloc (nonb->qpool);
}

#if 0
static void
p4est_dune_nquad_free (p4est_dune_nonb_t *nonb, p4est_quad_nonb_t *nquad)
{
  P4EST_ASSERT (nonb != NULL && nonb->qpool != NULL);
  P4EST_ASSERT (nquad != NULL);
  sc_mempool_free (nonb->qpool, nquad);
}
#endif

#ifndef P4_TO_P8
static const int    nonb_face_child[4][2] = {
  {0, 1},
  {2, 3},
  {0, 2},
  {1, 3},
};
#else
static const int    nonb_face_child[12][2] = {
  {0, 1},
  {2, 3},
  {4, 5},
  {6, 7},
  {0, 2},
  {1, 3},
  {4, 6},
  {5, 7},
  {0, 4},
  {1, 5},
  {2, 6},
  {3, 7},
};
#endif
static const int    nonb_face_number[3][2] = {
  {1, 0},
  {3, 2},
  {5, 4}
};

static void
p4est_dune_nonb_full (p4est_quad_nonb_t *nq, p4est_iter_face_side_t *fside)
{
  P4EST_ASSERT (nq != NULL);
  P4EST_ASSERT (fside != NULL);

  /* initialize a full-size face side */
  fside->is_hanging = 0;
  memset (&fside->is, 0, sizeof (fside->is));

  /* maybe there is one local quadrant */
  if (nq->squads.elem_count == 1) {
    P4EST_ASSERT (nq->ghosts.elem_count == 0);
    P4EST_ASSERT (nq->nvdesc == -1);

    /* this side is local */
    fside->is.full.is_ghost = 0;
    fside->is.full.quad = p4est_quadrant_array_index (&nq->squads, 0);
    fside->is.full.quadid = nq->quadid;
  }
  else {
    /* allow for ghost quadrants that should be there but aren't */
    P4EST_ASSERT (nq->squads.elem_count == 0);
    P4EST_ASSERT (nq->ghosts.elem_count <= 1);
    P4EST_ASSERT (nq->nvdesc != -1);

    /* this side is remote */
    fside->is.full.is_ghost = 1;
    if (nq->ghosts.elem_count == 0) {
      /* ghost not present when it should be by the mesh logic */
      P4EST_ASSERT (nq->nvdesc == 0);
      fside->is.full.quad = NULL;
      fside->is.full.quadid = -1;
    }
    else {
      /* proper full size ghost quadrant */
      P4EST_ASSERT (nq->nvdesc == -2);
      fside->is.full.quad = p4est_quadrant_array_index (&nq->ghosts, 0);
      fside->is.full.quadid = nq->ghostid;
    }
  }
}

static void
p4est_dune_nonb_hanging (p4est_quad_nonb_t *nq,
                         p4est_iter_face_side_t *fside, int ihang)
{
  P4EST_ASSERT (nq != NULL);
  P4EST_ASSERT (fside != NULL);

  /* initialize a hanging face side */
  fside->is_hanging = 1;
  memset (&fside->is, 0, sizeof (fside->is));

  /* maybe there is one local quadrant */
  if (nq->squads.elem_count == 1) {
    P4EST_ASSERT (nq->ghosts.elem_count == 0);
    P4EST_ASSERT (nq->nvdesc == -1);

    /* this side is local */
    fside->is.hanging.is_ghost[0] = 0;
    fside->is.hanging.quad[0] = p4est_quadrant_array_index (&nq->squads, 0);
    fside->is.hanging.quadid[0] = nq->quadid;
  }
  else {
    /* allow for ghost quadrants that should be there but aren't */
    P4EST_ASSERT (nq->squads.elem_count == 0);
    P4EST_ASSERT (nq->ghosts.elem_count <= 1);
    P4EST_ASSERT (nq->nvdesc != -1);

    /* this side is remote */
    fside->is.hanging.is_ghost[1] = 1;
    if (nq->ghosts.elem_count == 0) {
      /* ghost not present when it should be by the mesh logic */
      P4EST_ASSERT (nq->nvdesc == 0);
      fside->is.hanging.quad[0] = NULL;
      fside->is.hanging.quadid[0] = -1;
    }
    else {
      /* proper hanging ghost quadrant */
      P4EST_ASSERT (nq->nvdesc == -2);
      fside->is.hanging.quad[0] = p4est_quadrant_array_index (&nq->ghosts, 0);
      fside->is.hanging.quadid[0] = nq->ghostid;
    }
  }

  /* do not forget the hanging face number within parent face */
  P4EST_ASSERT (0 <= ihang && ihang < P4EST_HALF);
  fside->is.hanging.quadid[1] = ihang;
}

/* called with fully initialized nquad contexts for this interface */
static void
p4est_dune_nonb_face (p4est_dune_nonb_t *nonb,
                      p4est_quad_nonb_t *nquads[2], int ih)
{
  int                 i, j, k;
  int                 level[2];
  void              **lfound;
  p4est_iter_face_side_t *fside;
  p4est_quad_nonb_t  *nq, *fchildren[2][P4EST_HALF];
  p4est_quad_nonb_t  *cquads[2], nkey;

  /* for now, assume this is a face within a tree */
  P4EST_ASSERT (nonb != NULL);
  P4EST_ASSERT (nonb->finfo.orientation == 0);
  P4EST_ASSERT (nonb->finfo.sides.elem_count == 2);
  P4EST_ASSERT (nonb->iter_face != NULL);
  P4EST_ASSERT (nquads != NULL);
  P4EST_ASSERT (nquads[0] != NULL && nquads[1] != NULL);

  /* if there are no local quadrants on either side, we bail */
  if (nquads[0]->squads.elem_count == 0 && nquads[1]->squads.elem_count == 0) {
    return;
  }

  /* figure out if we are already in a hanging situation */
  for (k = 0; k < 2; ++k) {
    level[k] = nquads[k]->skey.level;
  }

  /* we will definitely execute the face callback or go into the recursion */
  for (k = 0; k < 2; ++k) {
    nq = nquads[k];
    fside = (p4est_iter_face_side_t *)
      sc_array_index (&nonb->finfo.sides, k);

    /* investigate the various possible cases */
    if (level[k] < level[!k]) {

      /* if this is a larger side, it must be whole and stays unchanged */
      P4EST_ASSERT (nq->nvdesc <= 0);
      P4EST_ASSERT (!fside->is_hanging);
    }
    else if (nq->nvdesc <= 0) {

      /* otherwise we must set new contents to the face side */
      P4EST_ASSERT (nq->squads.elem_count <= 1);
      P4EST_ASSERT (nq->ghosts.elem_count <= 1);

      /* same size is initalized as full */
      if (level[0] == level[1]) {
        p4est_dune_nonb_full (nq, fside);
      }
      else {
        P4EST_ASSERT (level[k] > level[!k]);

        /* initialize a hanging side */
        p4est_dune_nonb_hanging (nq, fside, ih);
      }
    }
    else {

      /* prepare what's needed to loop over the children below */
      for (j = 0; j < P4EST_HALF; ++j) {

        /* access the child quadrants touching the face */
        memset (&nkey, 0, sizeof (nkey));
        i = p4est_face_corners[fside->face][j];
        p4est_quadrant_child (&nq->skey, &nkey.skey, i);
        nkey.skey.p.which_tree = nq->skey.p.which_tree;
#if 0
        P4EST_LDEBUGF ("Looking up in tree %ld quadrant:\n",
                       (long) nkey.skey.p.which_tree);
        p4est_quadrant_print (SC_LP_ERROR, &nkey.skey);
#endif
        P4EST_EXECUTE_ASSERT_TRUE
          (sc_hash_lookup (nonb->qhash, &nkey, &lfound));
        fchildren[k][j] = (p4est_quad_nonb_t *) *lfound;
        P4EST_ASSERT (p4est_quadrant_child_id (&fchildren[k][j]->skey) == i);
      }
    }
  }

  /* if both sides are full size, we evaluate the face directly */
  if (nquads[0]->nvdesc <= 0 && nquads[1]->nvdesc <= 0) {
    P4EST_ASSERT (nquads[0]->squads.elem_count == 1 ||
                  nquads[1]->squads.elem_count == 1);

    /* execute the face callback with two unsplit sides */
    nonb->iter_face (&nonb->finfo, nonb->user_data);
  }
  else {
    /* go down the face recursion */
    for (j = 0; j < P4EST_HALF; ++j) {
      for (k = 0; k < 2; ++k) {
        cquads[k] = nquads[k]->nvdesc <= 0 ? nquads[k] : fchildren[k][j];
      }
      p4est_dune_nonb_face (nonb, cquads, j);
    }
  }
}

/* called with fully initialized nquad context for this volume */
static void
p4est_dune_nonb_volume (p4est_dune_nonb_t *nonb, p4est_quad_nonb_t *nquad)
{
  int                 i, j, k, m;
  size_t              nvz;
  size_t              oz, lz;
  size_t              goz, glz;
  void              **ifound;
  p4est_quad_nonb_t  *nchild;
  p4est_quad_nonb_t  *nchildren[P4EST_CHILDREN];
  p4est_quad_nonb_t  *nface[2];
  p4est_quadrant_t    children[P4EST_CHILDREN];
  p4est_quadrant_t   *tquad;
  p4est_iter_face_side_t *fside;

  /* verify preconditions */
  P4EST_ASSERT (nonb != NULL);
  P4EST_ASSERT (nquad != NULL);
  P4EST_ASSERT (p4est_quadrant_is_valid (&nquad->skey));
  P4EST_ASSERT (nonb->vinfo.treeid == nquad->skey.p.which_tree);

  /* handle special case of a full size leaf quadrant */
  if (nquad->squads.elem_count == 1) {
    tquad = p4est_quadrant_array_index (&nquad->squads, 0);

    /* if there is one real quadrant at this level, stop the recursion */
    if (tquad->level == nquad->skey.level) {
      /* mark full descendant as special case */
      nquad->nvdesc = -1;

      /* this is the place to execute the volume callback */
      if (nonb->iter_volume != NULL) {
        nonb->vinfo.quad = tquad;
        nonb->vinfo.quadid = nquad->quadid;
        nonb->iter_volume (&nonb->vinfo, nonb->user_data);
      }
      return;
    }
  }

  /* handle special case of a full size ghost quadrant */
  if (nquad->ghosts.elem_count == 1) {
    tquad = p4est_quadrant_array_index (&nquad->ghosts, 0);

    /* if there is one ghost quadrant at this level, stop the recursion */
    if (tquad->level == nquad->skey.level) {
      /* mark full descendant as special case */
      nquad->nvdesc = -2;
      return;
    }
  }

  /* now there is at least one quadrant below the current level */
  P4EST_ASSERT (nquad->split[P4EST_CHILDREN] == 0);
  p4est_split_array (&nquad->squads, nquad->skey.level, nquad->split);
  nvz = nquad->split[P4EST_CHILDREN];
  if (nonb->ghost_layer != NULL) {
    P4EST_ASSERT (nquad->gsplit[P4EST_CHILDREN] == 0);
    p4est_split_array (&nquad->ghosts, nquad->skey.level, nquad->gsplit);
    nvz += nquad->gsplit[P4EST_CHILDREN];
  }
  P4EST_ASSERT (nvz > 0);
  nquad->nvdesc = (p4est_locidx_t) nvz;

  /* loop through all children of current quadrant */
  p4est_quadrant_childrenv (&nquad->skey, children);
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    lz = nquad->split[i + 1] - (oz = nquad->split[i]);

    /* setup recursion quadrant object */
    nchild = nchildren[i] = p4est_dune_nquad_alloc (nonb);
    memset (nchild, 0, sizeof (*nchild));
    p4est_quadrant_copy (&children[i], &nchild->skey);
    nchild->skey.p.which_tree = nonb->vinfo.treeid;
    nchild->quadid = nquad->quadid + (p4est_locidx_t) oz;
    sc_array_init_view (&nchild->squads, &nquad->squads, oz, lz);

    /* setup ghost information similarly */
    glz = 0;
    if (nonb->ghost_layer != NULL) {
      glz = nquad->gsplit[i + 1] - (goz = nquad->gsplit[i]);
      nchild->ghostid = nquad->ghostid + (p4est_locidx_t) goz;
      sc_array_init_view (&nchild->ghosts, &nquad->ghosts, goz, glz);
    }

    /* hook the quadrant object into the hash table */
#if 0
    P4EST_LDEBUGF ("Inserting in tree %ld quadrant:\n",
                   (long) nchild->skey.p.which_tree);
    p4est_quadrant_print (SC_LP_ERROR, &nchild->skey);
#endif
    P4EST_EXECUTE_ASSERT_TRUE
      (sc_hash_insert_unique (nonb->qhash, nchild, &ifound));
    P4EST_ASSERT ((p4est_quad_nonb_t *) *ifound == nchild);

    /* we may have created a completely empty child context */
    if (lz + glz == 0) {
      continue;
    }

    /* we go into the general recursion algorithm */
    p4est_dune_nonb_volume (nonb, nchild);
  }

  /* call face recursion for all faces inside this volume */
  if (nonb->iter_face != NULL) {
    for (m = 0, i = 0; i < P4EST_DIM; ++i) {

      /* for a given dimension all faces face the same way */
      for (k = 0; k < 2; ++k) {
        fside = (p4est_iter_face_side_t *)
          sc_array_index (&nonb->finfo.sides, k);
        P4EST_ASSERT (nonb->vinfo.treeid == fside->treeid);
        fside->face = nonb_face_number[i][k];
        fside->is_hanging = 0;
      }

      /* go over all parallel faces in this direction */
      for (j = 0; j < P4EST_HALF; ++j, ++m) {
        for (k = 0; k < 2; ++k) {
          nface[k] = nchildren[nonb_face_child[m][k]];
        }
        p4est_dune_nonb_face (nonb, nface, j);
      }
    }
  }

  /* free temporary quadrant context */
  for (i = 0; i < P4EST_CHILDREN; ++i) {
    nchild = nchildren[i];

    /* Free memory here.  For tree faces, will have to defer this. */
#if 0
    P4EST_LDEBUGF ("Removing in tree %ld quadrant:\n",
                   (long) nchild->skey.p.which_tree);
    p4est_quadrant_print (SC_LP_ERROR, &nchild->skey);

    P4EST_EXECUTE_ASSERT_TRUE (sc_hash_remove (nonb->qhash, nchild, &rfound));
    P4EST_ASSERT ((p4est_quad_nonb_t *) rfound == nchild);
    p4est_dune_nquad_free (nonb, nchild);
#endif
  }
}

static void
p4est_dune_iterate_nonb (p4est_t * p4est, p4est_ghost_t * ghost_layer,
                         void *user_data, p4est_iter_volume_t iter_volume,
                         p4est_iter_face_t iter_face)
{
  int                 k;
  size_t              goz, glz;
  void              **ifound;
  p4est_topidx_t      tt;
  p4est_iter_face_side_t *fside;
  p4est_dune_nonb_t   snonb, *nonb = &snonb;
  p4est_quad_nonb_t  *nquad;

  P4EST_ASSERT (p4est != NULL);

  /* shortcuts */
  if (iter_volume == NULL && iter_face == NULL) {
    return;
  }

  /* setup context data */
  memset (nonb, 0, sizeof (*nonb));
  nonb->p4est = p4est;
  nonb->ghost_layer = ghost_layer;
  nonb->user_data = user_data;
  nonb->iter_volume = iter_volume;
  nonb->iter_face = iter_face;
  nonb->qhash = sc_hash_new (p4est_quad_nonb_hash,
                             p4est_quad_nonb_is_equal, nonb, NULL);
  nonb->qpool = sc_mempool_new (sizeof (p4est_quad_nonb_t));

  /* newly written hash cache */
  for (k = 0; k < P4EST_MAXLEVEL; ++k) {
    nonb->mru[k] = sc_hash_mru_new (p4est_quad_nonb_hash,
                                    p4est_quad_nonb_is_equal,
                                    NULL, NULL, 1 << 15);
  }

  /* prepare reusable volume context */
  nonb->vinfo.p4est = p4est;
  nonb->vinfo.ghost_layer = ghost_layer;
  nonb->vinfo.quad = NULL;
  nonb->vinfo.quadid = -1;
  nonb->vinfo.treeid = -1;

  /* prepare reusable face context */
  nonb->finfo.p4est = p4est;
  nonb->finfo.ghost_layer = ghost_layer;
  nonb->finfo.orientation = 0;
  nonb->finfo.tree_boundary = 0;
  sc_array_init_count
    (&nonb->finfo.sides, sizeof (p4est_iter_face_side_t), 2);

  /* loop over the local trees */
  for (tt = p4est->first_local_tree; tt <= p4est->last_local_tree; ++tt) {

    /* go into volume recursion with the tree root */
    nonb->tree = p4est_tree_array_index (p4est->trees, tt);
    nonb->vinfo.treeid = tt;
    for (k = 0; k < 2; ++k) {
      fside = (p4est_iter_face_side_t *)
        sc_array_index (&nonb->finfo.sides, k);
      fside->treeid = tt;
    }

    /* setup quadrant recursion state */
    nquad = p4est_dune_nquad_alloc (nonb);
    memset (nquad, 0, sizeof (*nquad));
    p4est_quadrant_root (&nquad->skey);
    nquad->skey.p.which_tree = tt;
    sc_array_init_view (&nquad->squads, &nonb->tree->quadrants,
                        0, nonb->tree->quadrants.elem_count);

    /* setup ghost information similarly */
    if (ghost_layer != NULL) {
      goz = nquad->ghostid = ghost_layer->tree_offsets[tt];
      glz = ghost_layer->tree_offsets[tt + 1] - goz;
      sc_array_init_view (&nquad->ghosts, &ghost_layer->ghosts, goz, glz);
    }

    /* hook the quadrant object into the hash table */
#if 0
    P4EST_LDEBUGF ("Inserting in tree %ld quadrant:\n",
                   (long) nquad->skey.p.which_tree);
    p4est_quadrant_print (SC_LP_ERROR, &nquad->skey);
#endif
    P4EST_EXECUTE_ASSERT_TRUE
      (sc_hash_insert_unique (nonb->qhash, nquad, &ifound));
    P4EST_ASSERT ((p4est_quad_nonb_t *) *ifound == nquad);

    /* we go into the general recursion algorithm */
    p4est_dune_nonb_volume (nonb, nquad);

    /* Free memory here.  For tree faces, will have to defer this. */
#if 0
    P4EST_LDEBUGF ("Removing in tree %ld quadrant:\n",
                   (long) nquad->skey.p.which_tree);
    p4est_quadrant_print (SC_LP_ERROR, &nquad->skey);
    P4EST_EXECUTE_ASSERT_TRUE (sc_hash_remove (nonb->qhash, nquad, &rfound));
    P4EST_ASSERT ((p4est_quad_nonb_t *) rfound == nquad);
    p4est_dune_nquad_free (nonb, nquad);
#endif
  }

  /* go into face recursion involving neighbor trees or boundary */
  if (nonb->iter_face != NULL) {

  }

  /* destroy hash cache */
  for (k = 0; k < P4EST_MAXLEVEL; ++k) {
    sc_hash_mru_destroy (nonb->mru[k]);
  }

  /* free context data */
  sc_hash_destroy (nonb->qhash);
  sc_mempool_destroy (nonb->qpool);
  sc_array_reset (&nonb->finfo.sides);
}

void
p4est_dune_iterate (p4est_t * p4est, p4est_ghost_t * ghost_layer,
                    void *user_data, p4est_iter_volume_t iter_volume,
                    p4est_iter_face_t iter_face)
{
  if (1) {
    p4est_dune_iterate_wrap (p4est, ghost_layer, user_data,
                             iter_volume, iter_face);
  }
  else {
    p4est_dune_iterate_nonb (p4est, ghost_layer, user_data,
                             iter_volume, iter_face);
  }
}
