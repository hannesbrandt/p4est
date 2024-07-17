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
#include <p4est_bits.h>
#include <p4est_dune.h>
#include <p4est_extended.h>
#include <p4est_vtk.h>
#else
#include <p8est_bits.h>
#include <p8est_dune.h>
#include <p8est_extended.h>
#include <p8est_vtk.h>
#endif /* P4_TO_P8 */

static int
refine_callback (p4est_t * p4est,
                 p4est_topidx_t which_tree, p4est_quadrant_t * quadrant)
{
  int                 maxlevel;
  int                 childid;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);

  /* do not refine beyond maximum level */
  maxlevel = *(int *) p4est->user_pointer;
  if (quadrant->level >= maxlevel) {
    return 0;
  }

  /* always refine a certain child for a tree */
  childid = p4est_quadrant_child_id (quadrant);
  if (childid == ((int) which_tree) % P4EST_CHILDREN) {
    return 1;
  }

  /* refine based on x coordinate */
  if (quadrant->x >= 9 * P4EST_QUADRANT_LEN (6) &&
      quadrant->x < 23 * P4EST_QUADRANT_LEN (6) &&
      (childid >= 1 && childid <= P4EST_HALF)) {
    return 1;
  }

  /* refine based on y coordinate */
  if (quadrant->y >= 5 * P4EST_QUADRANT_LEN (3) && childid >= 2) {
    return 1;
  }

#ifdef P4_TO_P8
  /* refine based on z coordinate */
  if (quadrant->z >= 6 * P4EST_QUADRANT_LEN (5) &&
      quadrant->z < 17 * P4EST_QUADRANT_LEN (5)) {
    return 1;
  }
#endif

  return 0;
}

static int
run_dune_interface (sc_MPI_Comm mpicomm, p4est_connectivity_t * conn,
                    int maxlevel, p4est_connect_type_t ctype)
{
  p4est_t            *p4est;
  p4est_ghost_t      *ghost;
  p4est_dune_numbers_t *dn;
  p4est_gloidx_t      gnum;
  int                 i;

  /* generate mesh with some arbitrary adaptive refinement */
  p4est = p4est_new_ext (mpicomm, conn, 0, 0, 1, 0, NULL, NULL);
  p4est->user_pointer = &maxlevel;

  /* run refinement loop */
  for (i = 0; i <= maxlevel; ++i) {
    gnum = p4est->global_num_quadrants;
    P4EST_GLOBAL_INFOF ("Into refinement iteration %d at %lld quadrants\n",
                        i, (long long) gnum);
    p4est_refine (p4est, 0, refine_callback, NULL);
    P4EST_ASSERT (gnum <= p4est->global_num_quadrants);
    if (gnum == p4est->global_num_quadrants) {
      break;
    }
    p4est_partition (p4est, 1, NULL);
  }
  P4EST_GLOBAL_INFOF ("After refinement iteration %d\n", i);

  /* wait with balancing to make the mesh more interesting */
  p4est_balance (p4est, P4EST_CONNECT_ALMOST, NULL);
  p4est_partition (p4est, 1, NULL);
  gnum = p4est->global_num_quadrants;
  P4EST_GLOBAL_INFOF ("Done after balance and partition at %lld quadrants\n",
                      (long long) gnum);

  /* output graphical representation */
  p4est_vtk_write_file (p4est, NULL, P4EST_STRING "_dune_interface");

  /* test various scenarios for DUNE node number export */
  for (i = 1; i < 2; ++i) {
    p4est_dune_numbers_params_t dparams, *pa = &dparams;
    p4est_dune_numbers_params_init (pa);
    pa->ctype = ctype;

    P4EST_GLOBAL_INFOF ("DUNE mesh interface iteration %d\n", i);

    /* we must provide a ghost layer */
    ghost = p4est_ghost_new (p4est, P4EST_CONNECT_FULL);

    /* generate node numbers for dune */
    dn = p4est_dune_numbers_new (p4est, ghost, pa);

    /* TO DO: do something with the DUNE node numbers */

    /* free memory in generated interface */
    p4est_dune_numbers_destroy (dn);

    /* deallocate temporary structure */
    p4est_ghost_destroy (ghost);
  }

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
  /* *INDENT-OFF* */
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
  /* *INDENT-ON* */
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
    if (run_dune_interface (mpicomm, conn, maxlevel, P4EST_CONNECT_FULL)) {
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
