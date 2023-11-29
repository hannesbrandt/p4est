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
#include <p4est_dune.h>
#include <p4est_extended.h>
#else
#include <p8est_dune.h>
#include <p8est_extended.h>
#endif /* P4_TO_P8 */

static int
run_dune_interface (sc_MPI_Comm mpicomm,
                    p4est_connectivity_t *conn, int maxlevel)
{
  p4est_t            *p4est;
  p4est_dune_numbers_t *dn;

  /* generate mesh with some arbitrary adaptive refinement */
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 1, 0, NULL, NULL);

  /* run refinement loop */

  /* generate node numbers for dune */
  dn = p4est_dune_numbers_new (p4est, NULL, NULL);
  p4est_dune_numbers_destroy (dn);

  /* deallocate generated mesh and return */
  p4est_destroy (p4est);
  return 0;
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  int                 progerr;
  int                 maxlevel;
  p4est_connectivity_t *conn;
  sc_MPI_Comm         mpicomm;
  const char         *usage;

  /* initialize MPI and p4est internals */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);
  mpicomm = sc_MPI_COMM_WORLD;

  sc_init (mpicomm, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  /* initialize global variables */
  conn = NULL;

  /* process command line arguments */
  progerr = 0;
  usage =
    "Arguments: <configuration> <level>\n"
    "   Configuration can be any of\n"
#ifndef P4_TO_P8
    "      unit|brick|three|moebius|star|periodic\n"
#else
    "      unit|brick|periodic|rotwrap|twocubes|rotcubes\n"
#endif
    "   Level is the maximum depth of refinement\n";
  if (!progerr && argc != 3) {
    P4EST_GLOBAL_LERROR ("Invalid command line argument count\n");
    progerr = 1;
  }
  if (!progerr) {
#ifndef P4_TO_P8
    if (!strcmp (argv[1], "unit")) {
      conn = p4est_connectivity_new_unitsquare ();
    }
    else if (!strcmp (argv[1], "brick")) {
      conn = p4est_connectivity_new_brick (2, 3, 0, 0);
    }
    else if (!strcmp (argv[1], "three")) {
      conn = p4est_connectivity_new_corner ();
    }
    else if (!strcmp (argv[1], "moebius")) {
      conn = p4est_connectivity_new_moebius ();
    }
    else if (!strcmp (argv[1], "star")) {
      conn = p4est_connectivity_new_star ();
    }
    else if (!strcmp (argv[1], "periodic")) {
      conn = p4est_connectivity_new_periodic ();
    }
#else
    if (!strcmp (argv[1], "unit")) {
      conn = p8est_connectivity_new_unitcube ();
    }
    else if (!strcmp (argv[1], "brick")) {
      conn = p8est_connectivity_new_brick (2, 3, 4, 0, 0, 0);
    }
    else if (!strcmp (argv[1], "periodic")) {
      conn = p8est_connectivity_new_periodic ();
    }
    else if (!strcmp (argv[1], "rotwrap")) {
      conn = p8est_connectivity_new_rotwrap ();
    }
    else if (!strcmp (argv[1], "twocubes")) {
      conn = p8est_connectivity_new_twocubes ();
    }
    else if (!strcmp (argv[1], "rotcubes")) {
      conn = p8est_connectivity_new_rotcubes ();
    }
#endif
    else {
      P4EST_GLOBAL_LERROR ("Invalid connectivity configuration\n");
      progerr = 1;
    }
  }
  if (!progerr) {
    maxlevel = atoi (argv[2]);
    if (maxlevel <= -1) {
      P4EST_GLOBAL_LERROR ("Invalid maximum level\n");
      progerr = 1;
    }
  }

  /* print usage message if error occured up to this point */
  if (progerr) {
    P4EST_GLOBAL_LERROR (usage);
  }

  /* run program proper */
  if (!progerr) {
    if (run_dune_interface (mpicomm, conn, maxlevel)) {
      P4EST_GLOBAL_LERROR ("Error running DUNE interface\n");
      progerr = 1;
    }
  }

  /* in the present version of this program the connectivity is not used */
  if (conn != NULL) {
    p4est_connectivity_destroy (conn);
  }

  /* clean up and exit */
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);
  return progerr ? EXIT_FAILURE : EXIT_SUCCESS;
}
