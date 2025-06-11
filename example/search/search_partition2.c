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

/** \file search_partition2.c
 *
 * This 2D example program searches random points in the partition of a 2D
 * brick forest and verifies the results by comparing them to the local search.
 */

#ifndef P4_TO_P8
#include <p4est_extended.h>
#include <p4est_search.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#endif

int
main (int argc, char **argv)
{
  int                 mpiret;
  p4est_connectivity_t *conn;
  p4est_t            *p4est;

  /* MPI initialization. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Package init. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* Create brick p4est. */
  conn = p4est_connectivity_new_brick (2, 2,
#ifdef P4_TO_P8
                                       2,
#endif
                                       0, 0
#ifdef P4_TO_P8
                                       , 0
#endif
    );
  p4est = p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, 3, 1, 0, NULL, NULL);

  /* Free memory. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Close MPI environment. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
