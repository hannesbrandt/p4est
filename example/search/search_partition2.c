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

#include <sc_options.h>
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
typedef struct search_partition_global
{
  /* p4est mesh */
  double              a[3], b[3], c[3]; /* refinement centers */
  int                 uniform_level;    /* level of initial uniform refinement */
  int                 max_level;        /* maximum level of adaptive refinement */
  p4est_t            *p4est;    /* the resulting p4est */

  /* query points */
  size_t              num_queries;      /* number of queries created on each process */
  int                 seed;     /* seed for random query creation */
  sc_array_t         *queries;  /* array of query points */

  /* search statistics */
  size_t              num_local_queries;        /* queries found in local search */
  sc_array_t         *global_nlq;       /* num_local_queries gathered globally */
}
search_partition_global_t;

typedef struct search_query
{
  double              xyz[3];   /* 3D coordinates */
  int                 is_local; /* set to 1, if found in local search */
  int                 rank;     /* rank assigned during partition search */
}
search_query_t;

static void
map_coordinates (p4est_qcoord_t x, p4est_qcoord_t y, p4est_qcoord_t z,
                 p4est_topidx_t which_tree, double *xyz)
{
  P4EST_ASSERT (which_tree < P4EST_CHILDREN);   /* assert we have a 2x2(x2) brick */
  xyz[0] = 0.5 * (COORDINATE_IROOTLEN * x + which_tree % 2);
  xyz[1] = 0.5 * (COORDINATE_IROOTLEN * y + (which_tree / 2) % 2);
#ifndef P4_TO_P8
  xyz[2] = 0.;
#else
  xyz[2] = 0.5 * (COORDINATE_IROOTLEN * z + which_tree / 4);
#endif
}

static int
refine_fn (p4est_t *p4est, p4est_topidx_t which_tree,
           p4est_quadrant_t *quadrant)
{
  p4est_qcoord_t      h2;
  double              xyz[3];
  double              dist, min_dist;

  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  search_partition_global_t *g =
    (search_partition_global_t *) p4est->user_pointer;
  P4EST_ASSERT (quadrant != NULL);

  /* get quadrant center reference coordinates in the unit square */
  P4EST_ASSERT (which_tree < P4EST_CHILDREN);   /* assert we have a 2x2(x2) brick */
  h2 = P4EST_QUADRANT_LEN (quadrant->level) >> 1;
  map_coordinates (quadrant->x + h2, quadrant->y + h2,
#ifndef P4_TO_P8
                   0,
#else
                   quadrant->z + h2,
#endif
                   which_tree, xyz);

  /* compute distance to point a */
  dist = (g->a[0] - xyz[0]) * (g->a[0] - xyz[0]) +
    (g->a[1] - xyz[1]) * (g->a[1] - xyz[1]) +
    (g->a[2] - xyz[2]) * (g->a[2] - xyz[2]);
  min_dist = sqrt (dist);

  /* compute distance to point b */
  dist = (g->b[0] - xyz[0]) * (g->b[0] - xyz[0]) +
    (g->b[1] - xyz[1]) * (g->b[1] - xyz[1]) +
    (g->b[2] - xyz[2]) * (g->b[2] - xyz[2]);
  min_dist = SC_MIN (min_dist, sqrt (dist));

  /* refine if quadrant center is close enough to either point a or point b */
  return (quadrant->level <
          g->max_level -
          floor (min_dist * (g->max_level - g->uniform_level) / 0.2));
}

static void
generate_queries (search_partition_global_t *g)
{
  size_t              iq, nqh;
  int                 id;
  search_query_t     *q;
  double              t;
  sc_array_t         *local_queries;

  /* generate local queries */
  local_queries =
    sc_array_new_count (sizeof (search_query_t), g->num_queries);
  sc_array_memset (local_queries, 0);
  nqh = local_queries->elem_count / 2;
  /* vary seeds between processes to get reproducible variety in points */
  srand (g->seed + 1000 * g->p4est->mpirank);
  for (iq = 0; iq < local_queries->elem_count; iq++) {
    q = (search_query_t *) sc_array_index (local_queries, iq);
    q->is_local = -1;
    q->rank = -1;
    for (id = 0; id < P4EST_DIM; id++) {
      q->xyz[id] = (double) rand () / RAND_MAX;
    }

    /* move point closer to g->b or g->c depending on iq and random t */
    t = pow ((double) rand () / RAND_MAX, 3);
    /* move the point to position sp->xyz * t + (1 - t) * {g->b,g->c} */
    if (iq < nqh) {
      q->xyz[0] = t * q->xyz[0] + (1 - t) * g->b[0];
      q->xyz[1] = t * q->xyz[1] + (1 - t) * g->b[1];
      q->xyz[2] = t * q->xyz[2] + (1 - t) * g->b[2];
    }
    else {
      q->xyz[0] = t * q->xyz[0] + (1 - t) * g->c[0];
      q->xyz[1] = t * q->xyz[1] + (1 - t) * g->c[1];
      q->xyz[2] = t * q->xyz[2] + (1 - t) * g->c[2];
    }
  }

  /* gather global queries on all processes */
  g->queries =
    sc_array_new_count (sizeof (search_query_t),
                        g->p4est->mpisize * g->num_queries);
  sc_array_memset (g->queries, 0);
  sc_MPI_Allgather (local_queries->array,
                    g->num_queries * sizeof (search_query_t), sc_MPI_BYTE,
                    g->queries->array,
                    g->num_queries * sizeof (search_query_t), sc_MPI_BYTE,
                    g->p4est->mpicomm);
  P4EST_GLOBAL_INFOF ("Created %ld global queries.\n",
                      g->queries->elem_count);

  /* cleanup */
  sc_array_destroy (local_queries);
}

static int
local_callback (p4est_t *p4est, p4est_topidx_t which_tree,
                p4est_quadrant_t *quadrant, p4est_locidx_t local_num,
                void *point)
{
  double              qxyz[3];
  double              qlen;
  double              tol;

  P4EST_ASSERT (point != NULL);
  search_query_t     *q = (search_query_t *) point;
  P4EST_ASSERT (p4est != NULL);
  P4EST_ASSERT (p4est->user_pointer != NULL);
  search_partition_global_t *g =
    (search_partition_global_t *) p4est->user_pointer;

  /* compute lower, left corners coords for quadrant bounds */
  map_coordinates (quadrant->x, quadrant->y,
#ifndef P4_TO_P8
                   0,
#else
                   quadrant->z,
#endif
                   which_tree, qxyz);
  qlen = 0.5 * P4EST_QUADRANT_LEN (quadrant->level) * COORDINATE_IROOTLEN;

  /* check if query is contained in quadrant */
  tol = 1e-14;
  if (q->xyz[0] < qxyz[0] - tol || q->xyz[0] > qxyz[0] + qlen + tol ||
      q->xyz[1] < qxyz[1] - tol || q->xyz[1] > qxyz[1] + qlen + tol ||
      q->xyz[2] < qxyz[2] - tol || q->xyz[2] > qxyz[2] + qlen + tol) {
    return 0;
  }

  if (local_num >= 0) {
    /* we are on a local leaf */
    q->is_local = 1;
    g->num_local_queries++;
  }

  return 1;
}

static void
search_local (search_partition_global_t *g)
{
  long long           lnq, *glnq, gnq;
  size_t              il;

  /* search queries locally */
  g->num_local_queries = 0;
  p4est_search_local (g->p4est, 0, NULL, local_callback, g->queries);
  P4EST_INFOF ("Queries found in local search = %ld\n", g->num_local_queries);

  /* allgather local num queries for future comparison with partition search */
  lnq = (long long) g->num_local_queries;
  glnq = P4EST_ALLOC (long long, g->p4est->mpisize);
  sc_MPI_Allgather (&lnq, 1, sc_MPI_LONG_LONG_INT, glnq, 1,
                    sc_MPI_LONG_LONG_INT, g->p4est->mpicomm);
  g->global_nlq = sc_array_new_count (sizeof (size_t), g->p4est->mpisize);
  gnq = 0;
  for (il = 0; il < g->global_nlq->elem_count; il++) {
    gnq += glnq[il];
    *(size_t *) sc_array_index (g->global_nlq, il) = glnq[il];
  }
  P4EST_GLOBAL_INFOF
    ("Queries found globally during local search = %lld (expected %ld)\n",
     gnq, g->queries->elem_count);
  P4EST_ASSERT (g->queries->elem_count <= (size_t) gnq);
  P4EST_FREE (glnq);
}

static void
run (search_partition_global_t *g)
{
  p4est_connectivity_t *conn;

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
  g->p4est =
    p4est_new_ext (sc_MPI_COMM_WORLD, conn, 0, g->uniform_level, 1, 0, NULL,
                   g);

  /* refine the forest adaptively around two points g->a and g->b */
  p4est_refine (g->p4est, 1, refine_fn, NULL);
  p4est_partition (g->p4est, 0, NULL);

  /* generate search queries */
  generate_queries (g);

  /* search queries in the local part of the mesh */
  search_local (g);

  /* output forest to vtk */
  p4est_vtk_write_file (g->p4est, NULL,
#ifndef P4_TO_P8
                        "search_partition2"
#else
                        "search_partition3"
#endif
    );

  /* Free memory. */
  sc_array_destroy (g->queries);
  sc_array_destroy (g->global_nlq);
  p4est_destroy (g->p4est);
  p4est_connectivity_destroy (conn);
}

int
main (int argc, char **argv)
{
  int                 mpiret;
  sc_options_t       *opt;
  int                 first_argc, ue;
  search_partition_global_t global, *g = &global;

  /* MPI initialization. */
  mpiret = sc_MPI_Init (&argc, &argv);
  SC_CHECK_MPI (mpiret);

  /* Package init. */
  sc_init (sc_MPI_COMM_WORLD, 1, 1, NULL, SC_LP_DEFAULT);
  p4est_init (NULL, SC_LP_DEFAULT);

  opt = sc_options_new (argv[0]);
  sc_options_add_int (opt, 'u', "uniform_level", &g->uniform_level, 3,
                      "Level of uniform refinement");
  sc_options_add_int (opt, 'm', "max_level", &g->max_level, 7,
                      "Level of maximum refinement");
  sc_options_add_size_t (opt, 'q', "num_queries", &g->num_queries, 100,
                         "Number of queries created per process");
  sc_options_add_int (opt, 's', "seed", &g->seed, 0,
                      "Seed for random queries");

  /* proceed in run-once loop for clean abort */
  ue = 0;
  do {
    first_argc = sc_options_parse (p4est_package_id, SC_LP_DEFAULT,
                                   opt, argc, argv);
    if (first_argc < 0) {
      P4EST_GLOBAL_LERROR ("Invalid option format.\n");
      ue = 1;
      break;
    }
    sc_options_print_summary (p4est_package_id, SC_LP_ESSENTIAL, opt);

    /* check options for consistency */
    if (g->uniform_level < 0 || g->uniform_level > P4EST_QMAXLEVEL) {
      P4EST_GLOBAL_LERRORF ("Uniform level out of bounds 0..%d\n",
                            P4EST_QMAXLEVEL);
      ue = 1;
    }
    if (g->max_level < 0 || g->max_level > P4EST_QMAXLEVEL) {
      P4EST_GLOBAL_LERRORF ("Maximum level out of bounds 0..%d\n",
                            P4EST_QMAXLEVEL);
      ue = 1;
    }
    if (g->seed < 0) {
      P4EST_GLOBAL_LERROR ("Seed has to be non-negative.\n");
      ue = 1;
    }
    if (ue) {
      break;
    }

    /* define centers of refinement and point creation */
    g->a[0] = 0.2;
    g->a[1] = 0.4;
    g->a[2] = 0.4;
    g->b[0] = 0.7;
    g->b[1] = 0.55;
    g->b[2] = 0.55;
    g->c[0] = 0.3;
    g->c[1] = 0.6;
    g->c[2] = 0.8;
#ifndef P4_TO_P8
    g->a[2] = g->b[2] = g->c[2] = 0.;   /* reset z-coordinate to 0 */
#endif

    /* run example */
    run (g);
  }
  while (0);
  if (ue) {
    sc_options_print_usage (p4est_package_id, SC_LP_ERROR, opt, NULL);
  }

  /* Close MPI environment. */
  sc_options_destroy (opt);
  sc_finalize ();
  mpiret = sc_MPI_Finalize ();
  SC_CHECK_MPI (mpiret);

  return EXIT_SUCCESS;
}
