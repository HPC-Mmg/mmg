/* =============================================================================
**  This file is part of the mmg software package for the tetrahedral
**  mesh modification.
**  Copyright (c) Bx INP/CNRS/Inria/UBordeaux/UPMC, 2004-
**
**  mmg is free software: you can redistribute it and/or modify it
**  under the terms of the GNU Lesser General Public License as published
**  by the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  mmg is distributed in the hope that it will be useful, but WITHOUT
**  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
**  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
**  License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License and of the GNU General Public License along with mmg (in
**  files COPYING.LESSER and COPYING). If not, see
**  <http://www.gnu.org/licenses/>. Please read their terms carefully and
**  use this copy of the mmg distribution only if you accept them.
** =============================================================================
*/

/**
 * \file octree4immersedbdy_3d.c
 * \brief Creation and refinment of an octree mesh to capture a surface.
 * \author Gustavo Borgez Valim (Enseirb-Matmeca)
 * \author Laurent Le (Enseirb-Matmeca)
 * \author Damien Sans (Enseirb-Matmeca)
 * \author Algiane Froehly (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Creation and refinment of an octree mesh to capture an immersed boundary.
 *
 */

#include "mmg3d.h"
#include "../mmgs/mmgs.h"

#define MMG3D_IMMERSRAT 2 // ratio between the octree and the immersed surface

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Create the initial octree cell that will embed the provided surface (stored
 * in mesh->tria, mesh->points...). We suppose here that the octree root will have a
 * size of MMG3D_IMMERSRAT \f$ \times \f$ the surface mesh bounding box and that
 * the immersed surface is centered inside the octree root.
 *
 * \remark the surface mesh must be scaled so the bounding box info have been
 * computed.
 *
 */
static inline
int MMG3D_create_bbOctree(MMG5_pMesh mesh, MMG5_pSol sol) {

  /* if mesh->info.delta==0 it means that the BB has not been computed */
  assert ( mesh->info.delta > 0.  );

  printf("TO IMPLEMENT\n");
  return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Refine the octree in order to have only one point per leaf.
 *
 */
static inline
int MMG3D_refine_octreeOnPoints(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf("TO IMPLEMENT\n");
  return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Refine the octree in order to have:
 *  - if the octree leaf contains a point: only triangles that contains this point
 *  - if the octree leaf doesn't contains any point: only 1 triangle that
 *     intersect the leaf
 *
 */
static inline
int MMG3D_refine_octreeOnTria(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf("TO IMPLEMENT\n");
  return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Balance the octree (2:1) to have at max 1 level of refinement between 2
 * adjacent cells
 *
 */
static inline
int MMG3D_balance_octree(MMG5_pMesh mesh, MMG5_pSol sol) {

  printf("TO IMPLEMENT\n");
  return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure (possibly empty).
 *
 * \return 1 if success, 0 if fail.
 *
 * Creation of an adapted octree mesh to capture an immersed boundary.
 *
 */
int MMG3D_octree_for_immersedBdy(MMG5_pMesh mesh, MMG5_pSol sol) {

  /**--- stage 1: Immersed surface analysis */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** IMMERSED SURFACE ANALYSIS\n");

  /* Scale the input surface mesh between [0;1] (along the largest direction)
   * and compute the immersed mesh bounding box */
  if ( !MMG5_scaleMesh( mesh, NULL,sol) ) {
    fprintf(stderr,"\n  ## Mesh scaling problem. Exit program.\n");
    return 0;
  }

  /* Compute adjacency relationship */
  if ( !MMGS_hashTria( mesh) ) {
    fprintf(stderr,"\n  ## Hashing problem. Exit program.\n");
    return 0;
  }

  /**--- stage 2: Octree initialization */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** OCTREE INITIALIZATION\n");

  /* Creation of the octree root  */
  if ( !MMG3D_create_bbOctree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree root creation problem. Exit program.\n");
    return 0;
  }
#ifndef NDEBUG
  // To remove when the octree creation will be ok
  if ( !MMG3D_saveVTKOctree(mesh,sol,"bbOctree.vtk") ) {
    fprintf(stderr,"\n  ## Warning: unable to save the initial octree\n");
  }
#endif

  /**--- stage 3: Octree refinement */
  if ( abs(mesh->info.imprim) > 3 )
    fprintf(stdout,"\n  ** OCTREE REFINMENT\n");

  /* Octree refinement in order to have 1 node per leaf */
  if ( !MMG3D_refine_octreeOnPoints(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree refinement (over points) problem."
            " Exit program.\n");
    return 0;
  }

  /* Octree refinement in order to have 1 boundary entity (triangle) per leaf */
  if ( !MMG3D_refine_octreeOnTria(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree refinement (over triangles) problem."
            " Exit program.\n");
    return 0;
  }

  /* Octree refinement in order to have 1 boundary entity (triangle) per leaf */
  if ( !MMG3D_refine_octreeOnTria(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree refinement (over triangles) problem."
            " Exit program.\n");
    return 0;
  }

  /* 2:1 balancing of the octree  */
  if ( !MMG3D_balance_octree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree refinement (over triangles) problem."
            " Exit program.\n");
    return 0;
  }

  /* Unscale the surface mesh and the octree */
  if ( !MMG5_unscaleMesh( mesh, NULL,sol) ) {
    fprintf(stderr,"\n  ## Mesh scaling problem. Exit program.\n");
    return 0;
  }

#ifndef NDEBUG
  // To remove when the octree creation will be ok
  if ( !MMG3D_saveVTKOctree(mesh,sol,"finalOctree.vtk") ) {
    fprintf(stderr,"\n  ## Warning: unable to save the initial octree\n");
  }
#endif

  /* Memory release */
  MMG3D_delete_octree (mesh);

  return 1;
}
