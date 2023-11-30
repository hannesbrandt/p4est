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

/** \file p8est_dune.h
 * Populate element corner, edge, and face numbers to interface to DUNE.
 * Depending on the invocation, they are locally or globally unique.
 * The same number may be used for some corner, some edge and some face.
 * Uniqueness holds among corners, separately among edges and among faces.
 */

#ifndef P8EST_DUNE_H
#define P8EST_DUNE_H

#include <p8est_ghost.h>

SC_EXTERN_C_BEGIN;

/** Parameters to the constructor of the dune number tables. */
typedef struct p8est_dune_numbers_params
{
  /** Determine which codimensions are numbered. */
  p8est_connect_type_t ctype;
}
p8est_dune_numbers_params_t;

/** Element corner, edge and face numbers to interface to the DUNE library.
 *
 * The element corner, edge and face arrays hold global indices
 * that are unique within the partition among all faces, and likewise
 * unique within the partition among all corners, and among all edges.
 *
 * If the \ref p8est_ghost_t member to \ref p8est_dune_numbers_new has
 * been passed as NULL, then the numbers will not be globally synced.
 *
 * The global indices are of type \ref p4est_gloidx_t.
 */
typedef struct p8est_dune_numbers
{
  /** Parameters the structure was built with copied by value. */
  p8est_dune_numbers_params_t params;

  /** For each element corner a partition-independent global index. */
  sc_array_t         *element_corners;

  /** For each element edge a partition-independent global index. */
  sc_array_t         *element_edges;

  /** For each element face a partition-independent global index. */
  sc_array_t         *element_faces;
}
p8est_dune_numbers_t;

/** Set default parameters to pass to \ref p8est_dune_numbers_new.
 * \param [out] params  Pointer must not be NULL.
 *                      The structure is filled with default values.
 */
void                p8est_dune_numbers_params_default
  (p8est_dune_numbers_params_t * params);

/** Create lookup tables for unique corners, edges and faces.
 * Hanging corners, edges and faces are included.
 * All element corner numbers are partition-independent und mutually distinct.
 * All element edge numbers are partition-independent und mutually distinct.
 * All element face numbers are partition-independent und mutually distinct.
 * \param [in] p4est    Required input parameter is not modified.
 * \param [in] ghost    The ghost layer must match the \a p4est passed.
 *                      It may also be NULL, in which case the numbers are
 *                      local to the process and not globally meaningful.
 * \param [in] params   Further parameters to control the mode of operation.
 *                      When passing NULL, the behavior is identical to using
 *                      defaults by \ref p8est_dune_numbers_params_default.
 * \return              A fully initialized dune node numbering.
 *                      Deallocate with \ref p8est_dune_numbers_destroy.
 */
p8est_dune_numbers_t *p8est_dune_numbers_new (p8est_t * p4est,
                                              p8est_ghost_t * ghost,
                                              const
                                              p8est_dune_numbers_params_t *
                                              params);

/** Destroy the dune element number tables previously generated.
 * \param [in] dn       Valid dune number tables from \ref
 *                      p8est_dune_numbers_new are deallocated.
 */
void                p8est_dune_numbers_destroy (p8est_dune_numbers_t * dn);

SC_EXTERN_C_END;

#endif /* P8EST_DUNE_H */
