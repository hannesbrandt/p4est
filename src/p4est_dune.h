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

/** \file p4est_dune.h
 *
 * Populate element corner and face number tables to interface to DUNE.
 * The same number may be used for some corner and some face.
 * Uniqueness holds among corners, and separately among faces.
 *
 * \ingroup p4est
 */

#ifndef P4EST_DUNE_H
#define P4EST_DUNE_H

#include <p4est_ghost.h>

SC_EXTERN_C_BEGIN;

/** Parameters to the constructor of the dune number tables. */
typedef struct p4est_dune_numbers_params
{
  /** Determine which codimensions are numbered. */
  p4est_connect_type_t ctype;
}
p4est_dune_numbers_params_t;

/** Element corner and face numbers to interface to the DUNE library.
 *
 * The element corners and faces arrays hold local indices that are unique
 * within this process among all faces, and likewise unique within this
 * process among all corners.  The same number may occur in both the
 * face and corner list, but it refers to a different mesh object.
 *
 * The indices are numbered starting from zero.
 *
 * The array entries are of type \ref p4est_locidx_t.
 */
typedef struct p4est_dune_numbers
{
  /** Parameters the numbers are built with are copied by value. */
  p4est_dune_numbers_params_t params;

  /** Highest number (exclusive) among all array entries. */
  p4est_locidx_t      num_local_numbers;

  /** For each element corner a process-local index. */
  sc_array_t         *element_corners;

  /** For each element face a process-local index. */
  sc_array_t         *element_faces;
}
p4est_dune_numbers_t;

/** Set default parameters to pass to \ref p4est_dune_numbers_new.
 * \param [out] params  Pointer must not be NULL.
 *                      The structure is filled with default values.
 */
void                p4est_dune_numbers_params_default
  (p4est_dune_numbers_params_t * params);

/** Create lookup tables for unique corners and faces.
 * Hanging corners and faces are included as independent entities.
 * All numbers are process-local.
 *
 * \param [in] p4est    Required input parameter is not modified.  It must
 *                      be balanced at least to \ref P4EST_CONNECT_ALMOST.
 * \param [in] ghost    The ghost layer must have been generated with
 *                      \ref p4est_ghost_new using the same \c p4est and
 *                      the parameter \ref P4EST_CONNECT_FULL.
 * \param [in] params   Further parameters to control the mode of operation.
 *                      When passing NULL, the behavior is identical to using
 *                      defaults by \ref p4est_dune_numbers_params_default.
 * \return              A fully initialized dune node numbering.
 *                      Deallocate with \ref p4est_dune_numbers_destroy.
 */
p4est_dune_numbers_t *p4est_dune_numbers_new (p4est_t * p4est,
                                              p4est_ghost_t * ghost,
                                              const
                                              p4est_dune_numbers_params_t *
                                              params);

/** Destroy the dune element number tables previously generated.
 * \param [in] dn       Valid dune number tables from \ref
 *                      p4est_dune_numbers_new are deallocated.
 */
void                p4est_dune_numbers_destroy (p4est_dune_numbers_t * dn);

SC_EXTERN_C_END;

#endif /* P4EST_DUNE_H */
