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
#include <p4est_dune.h>
#include <p4est_lnodes.h>
#else
#include <p8est_algorithms.h>
#include <p8est_dune.h>
#include <p8est_lnodes.h>
#endif

/* map element corner to element node index */
static const int    corner_enode[P4EST_CHILDREN] = {
  0, 4, 20, 24
#ifdef P4_TO_P8
    , 100, 104, 120, 124
#endif
};

/* set default parameters */
void
p4est_dune_numbers_params_default (p4est_dune_numbers_params_t * params)
{
  P4EST_ASSERT (params != NULL);

  /* set default parameters to construct a dune numbers structure */
  params->ctype = P4EST_CONNECT_FULL;
}

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

/* Populate number arrays for real.  This function does all the work. */
static void
generate_numbers (p4est_dune_numbers_t * dn, p4est_lnodes_t * ln)
{
  p4est_lnodes_code_t fc;
  p4est_locidx_t      lne, e;
  p4est_locidx_t      lni;
  p4est_locidx_t     *eno;
  p4est_gloidx_t     *geco;
  int                 c;
  int                 do_c, do_f;
#ifdef P4_TO_P8
  int                 do_e;
#endif

  /* basic verification of input parameters */
  P4EST_ASSERT (dn != NULL);
  P4EST_ASSERT (ln != NULL);
  P4EST_ASSERT (ln->degree == 4);

  /* examine which codimensions are handled */
  do_c = dn->element_corners != NULL;
#ifdef P4_TO_P8
  do_e = 0;
#endif
  do_f = 0;
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
  eno = ln->element_nodes;
  for (e = 0; e < lne; ++e) {
    /* retrieve pointers to place output */
    if (do_c) {
      geco = (p4est_gloidx_t *)
        sc_array_index (dn->element_corners, e * P4EST_CHILDREN);
    }

    /* check for hanging status of the element */
    if (!(fc = ln->face_code[e])) {
      /* the element is not the small neighbor to any other */
      if (do_c) {
        for (c = 0; c < P4EST_CHILDREN; ++c) {
          /* retrieve local node number */
          lni = eno[corner_enode[c]];
          geco[c] = lni_to_gni (ln, lni);
        }
      }
      if (do_f) {
        /* face nodes not yet implemented */
        SC_ABORT_NOT_REACHED ();
      }
    }
    else {
      /* hanging elements are not yet implemented */
      SC_ABORT_NOT_REACHED ();
    }
    eno += P4EST_CHILDREN;
  }
}

/* create lookup tables for unique corners and faces */
p4est_dune_numbers_t *
p4est_dune_numbers_new (p4est_t * p4est, p4est_ghost_t * ghost,
                        const p4est_dune_numbers_params_t * params)
{
  p4est_dune_numbers_t *dn;
  p4est_dune_numbers_params_t *pa;
  p4est_ghost_t      *gh;
  p4est_lnodes_t     *ln;
  p4est_locidx_t      lne;

  /* formally verify arguments */
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est_is_valid (p4est));

  /* allocate numbers structure to populate */
  dn = P4EST_ALLOC_ZERO (p4est_dune_numbers_t, 1);
  pa = &dn->params;

  /* establish input parameters and remember them */
  if (params == NULL) {
    p4est_dune_numbers_params_default (pa);
  }
  else {
    *pa = *params;
    params = NULL;
  }

  /* generate temporary local ghost (if needed) and node numbering */
  ln = p4est_lnodes_new
    (p4est, ghost == NULL ?
     (gh = p4est_ghost_new_local (p4est, pa->ctype)) : (gh = NULL, ghost), 4);
  if (gh != NULL) {
    /* this was a local lnodes without parallel neighbors */
    p4est_ghost_destroy (gh);
  }

  /* construct corner and face numbers and return them */
  lne = p4est->local_num_quadrants;
  if (pa->ctype >= P4EST_CONNECT_CORNER) {
    dn->element_corners =
      sc_array_new_count (sizeof (p4est_gloidx_t), lne * P4EST_CHILDREN);
  }
  if (pa->ctype >= P4EST_CONNECT_FACE) {
    dn->element_faces =
      sc_array_new_count (sizeof (p4est_gloidx_t), lne * P4EST_FACES);
  }
  generate_numbers (dn, ln);

  /* clean internal state and return */
  p4est_lnodes_destroy (ln);
  return dn;
}

void
p4est_dune_numbers_destroy (p4est_dune_numbers_t * dn)
{
  P4EST_ASSERT (dn != NULL);

  if (dn->element_corners != NULL) {
    sc_array_destroy (dn->element_corners);
  }
  if (dn->element_faces != NULL) {
    sc_array_destroy (dn->element_faces);
  }
  P4EST_FREE (dn);
}
