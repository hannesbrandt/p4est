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
#include <p4est_vtk.h>
#else
#include <p8est_extended.h>
#include <p8est_search.h>
#include <p8est_vtk.h>
#endif

#define COORDINATE_IROOTLEN (1. / P4EST_ROOT_LEN)

/* define centers of refinement and point creation */
#ifndef P4_TO_P8
const double        point_a[3] = { 0.2, 0.4, 0. };
const double        point_b[3] = { 0.7, 0.5, 0. };
const double        point_c[3] = { 0.3, 0.6, 0. };
#else
const double        point_a[3] = { 0.2, 0.4, 0.4 };
const double        point_b[3] = { 0.7, 0.5, 0.5 };
const double        point_c[3] = { 0.3, 0.6, 0.8 };
#endif

static int
refine_fn (p4est_t *p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t *quadrant)
{
  p4est_qcoord_t      h2;
  double              xyz[3];
  double              dist, min_dist;

  /* get quadrant center reference coordinates in the unit square */
  P4EST_ASSERT (which_tree < P4EST_CHILDREN);   /* assert we have a 2x2(x2) brick */
  h2 = P4EST_QUADRANT_LEN (quadrant->level) >> 1;
  xyz[0] = 0.5 * (COORDINATE_IROOTLEN * (quadrant->x + h2) + which_tree % 2);
  xyz[1] =
    0.5 * (COORDINATE_IROOTLEN * (quadrant->y + h2) + (which_tree / 2) % 2);
#ifndef P4_TO_P8
  xyz[2] = 0.;
#else
  xyz[2] = 0.5 * (COORDINATE_IROOTLEN * (quadrant->z + h2) + which_tree / 4);
#endif

  /* compute distance to point a */
  dist = (point_a[0] - xyz[0]) * (point_a[0] - xyz[0]) +
    (point_a[1] - xyz[1]) * (point_a[1] - xyz[1]) +
    (point_a[2] - xyz[2]) * (point_a[2] - xyz[2]);
  min_dist = sqrt (dist);

  /* compute distance to point b */
  dist = (point_b[0] - xyz[0]) * (point_b[0] - xyz[0]) +
    (point_b[1] - xyz[1]) * (point_b[1] - xyz[1]) +
    (point_b[2] - xyz[2]) * (point_b[2] - xyz[2]);
  min_dist = SC_MIN (min_dist, sqrt (dist));

  /* refine if quadrant center is close enough to either point a or point b */
  return (quadrant->level < 7 - floor (min_dist / 0.05));
}

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
  p4est_refine (p4est, 1, refine_fn, NULL);

  /* output forest to vtk */
  p4est_vtk_write_file (p4est, NULL,
#ifndef P4_TO_P8
                        "search_partition2"
#else
                        "search_partition3"
#endif
    );

  /* Free memory. */
  p4est_destroy (p4est);
  p4est_connectivity_destroy (conn);

  /* Close MPI environment. */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
