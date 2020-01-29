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

#define MMG3D_IMMERSRAT 0.5 // ratio between the immersed surface and the octree

/**
 * \param mesh pointer toward a mesh structure.
 * \param sol pointer toward a solution structure.
 *
 * \return 1 if success, 0 if fail.
 *
 * Create the initial octree cell that will embed the provided surface (stored
 * in mesh->tria, mesh->points...). We suppose here that the octree root will
 * represent the [0;1]^3 cell.
 *
 * \remark the surface mesh must be scaled, thus, the bounding box info is
 * already computed.
 *
 */
static inline
int MMG3D_create_bbOctree(MMG5_pMesh mesh, MMG5_pSol sol) {

  /* if mesh->info.delta==0 it means that the BB has not been computed */
  assert ( mesh->info.delta > 0.  );

  /* Initialization of the first octree cell. */
  if ( !MMG3D_init_MIBOctree(mesh,&mesh->iboctree) )
    return 0;


  return 1;
}


/**
 * \param a the Z-value of a mesh point
 * \param b the Z-value of another mesh point
 * \param mesh pointer toward a mesh structure.
 * \param q pointer toward the MOctree cell
 * \param root pointer toward the octree root
 *
 * \return 1 if success, 0 if fail.
 *
 * Recursive algorithm used to refine the octree
 *
 */
static inline
int MMG3D_recursive_refine(int64_t a, int64_t b, MMG5_pMesh mesh, MMG5_MIBOctree_s *q, MMG5_pMIBOctree *root) {
  int i = 0;
  int64_t comp_value = 0;
  comp_value = (int64_t) 7 << (3*((*root)->depth_max - q->depth - 1));   // ie 111 << 2^{n}
  // compare the trio 000 at the same depth for a and b
  if ((q->depth < (*root)->depth_max - 1 ) && ((a & comp_value) == ( b & comp_value))  )
  {
    int j = 0;
    MMG5_MIBOctree_s* in_son = mesh->iboctree->root + q->sons[0];

    // Find the son containing the points
    while ((j < 7) && ((in_son->Z_coord & comp_value) != ( a & comp_value)))
    {
      j += 1;
      in_son = mesh->iboctree->root + q->sons[j];
    }
    // If the does not already have any sons, split the cell and recover the maximal reached depth
    if (in_son->is_filled == 0){
      MMG3D_split_MIBOctree_s ( mesh,in_son,&mesh->iboctree );
      in_son->is_filled = 1 ;
      if (mesh->iboctree->max_depth_reached < in_son->depth)
        mesh->iboctree->max_depth_reached = in_son->depth + 1;
    }
    MMG3D_recursive_refine( a, b, mesh, in_son, &mesh->iboctree);
  }
  return 1;
}


/**

 * \param x x-coordinates of a mesh point
 * \param y y-coordinates of a mesh point
 * \param z z-coordinates of a mesh point
 *
 * \return Z a 64-bit integer.
 *
 * Compute the Z-value of a mesh point given its coordinates
 *
 */
static inline
int64_t coords_into_Z_coord(int x, int y, int z) {
  if (x > 4194303 || y > 4194303 || z > 4194303)
    printf("\n %s:%d:\n%s: Warning ! Coordinates may be too big to be converted into Z-order\n",__FILE__,__LINE__,__func__);

  int i = 0, j = 0;

  int64_t Z = 0;
  for (i = 0 ; i < 22 ; i++)
  {
    j = 2*i;
    Z += (int64_t)(x & (1<<i)) << j;
    j++;
    Z += (int64_t)(y & (1<<i)) << j;
    j++;
    Z += (int64_t)(z & (1<<i)) << j;
  }
  return Z;
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

  int i = 0;
  int j = 0;
  int64_t *points_Z = (int64_t*) calloc(mesh->np, sizeof(int64_t));
  int *scaled_coord = (int*) calloc(3, sizeof(int));

  // Point scaling
  for (i = 1 ; i <= mesh->np ; i++)
  {
    for (j = 0; j < 3 ; j++)
    {
      scaled_coord[j] = (int) (mesh->point[i].c[j]*(1 << mesh->iboctree->depth_max)); // multiply by 2^{depth_max}
    }
    points_Z[i-1] = coords_into_Z_coord(scaled_coord[0], scaled_coord[1], scaled_coord[2]);
  }

  // Split the root
  MMG3D_split_MIBOctree_s ( mesh,&mesh->iboctree->root[0],&mesh->iboctree );

  // Call the recursive function rolling over every mesh point
  for (i = 0 ; i < mesh->np - 1 ; i++)
  {
    for (j = i+1 ; j < mesh->np ; j++)
    {
      MMG3D_recursive_refine( points_Z[i], points_Z[j], mesh, &mesh->iboctree->root[0], &mesh->iboctree);
    }

  }
  free(points_Z);
  free(scaled_coord);

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

  printf("\n %s:%d:\n%s: FUNCTION TO IMPLEMENT\n",__FILE__,__LINE__,__func__);
  return 0;

  return 1;
}

/**
 * \param mesh pointer toward a mesh structure.
 * \param root pointer toward the octree root
 * \param Z a Z value
 *
 * \return in_son pointer towards a cell.
 *
 * Find the deepest cell corresponding to the given Z-value
 *
 */
static inline
MMG5_MIBOctree_s* MMG3D_find_octree(MMG5_pMesh mesh, MMG5_pMIBOctree *root, int64_t Z) {

  int i = 0;
  int64_t comp_value = 0;
  MMG5_MIBOctree_s *in_son = &mesh->iboctree->root[0] ;

  while ( i < (*root)->depth_max && in_son->sons[0] != 0 ){
    comp_value = (int64_t) 7 << (3*((*root)->depth_max - i - 1));   // ie 111 << 2^{n}
    comp_value = (Z & comp_value) >> (3*((*root)->depth_max - i - 1));
    in_son = mesh->iboctree->root + in_son->sons[comp_value];
    i += 1;
  }

  return (in_son);
}



/**
 * \param root pointer toward the octree root
 * \param depth the depth value of the cell we look at
 * \param direction integer entailing the axe direction (0 for x, 1 for y, 2 for z)
 *
 * \return Z_n a Z-value corresponding to a neighbouring cell
 *
 * Find the Z value of the adjacent cell in the x/y/z direction
 *
 */
 int64_t MMG3D_external_neighbours(MMG5_pMIBOctree *root, int64_t Z, int depth, int direction){

   int i;
   int64_t comp_value;
   int64_t a = 0;
   int64_t b = 0;
   int deplacement;
   int64_t Z_n = (int64_t) 0;

     deplacement = 3*((*root)->depth_max - depth);
     comp_value = (int64_t) 1 << (deplacement + direction);
     i = depth;
     while ( i > 0 &&  a == b){
       a = b;
       if ((Z & comp_value) == comp_value){
         b = 1;
       }else if((Z & comp_value)!= comp_value){
         b = 0;
       }else{
         printf("ERROR : Comparative value non valid ! l.%d\n",__LINE__);
       }
       if ( i == depth)
         a = b;

       if (direction == 0){
         Z_n += (b^1) << deplacement;
         Z_n += (Z & ((int64_t)3 << (deplacement + 1)));

       }else if (direction == 1){
         Z_n += (Z & ((int64_t)1 << deplacement));
         Z_n += (b^1) << (deplacement+1);
         Z_n += (Z & ((int64_t)1 << (deplacement+2)));

       }else if (direction == 2){
         Z_n += (Z & ((int64_t)3 << (deplacement)));
         Z_n += (b^1) << (deplacement+2);
       }else{
         printf(" WARNING : no orientation corresponding to this choice of direction : l.%d\n",__LINE__);
       }

       comp_value = comp_value << 3;
       deplacement += 3;
       i -= 1;
     }
     if (i == 0 && a == b){
       Z_n = -1;
     }else{
       while ( i > 0 ){
         Z_n += (Z & ((int64_t)7 << (deplacement)));
         deplacement += 3;
         i -= 1;
       }
     }

   return (Z_n);
 }


/**
 * \param mesh pointer toward a mesh structure.
 * \param q pointer toward the MOctree cell
 * \param root pointer toward the octree root
 * \param Z 64-bit integer guiding the refinement
 * \param deepest the depth required for the original cell to be balanced
 *
 * \return 1 if success, 0 if fail.
 *
 * Recursive algorithm splitting cells until the original cell is balanced
 *
 */
static inline
int MMG3D_split_in_balance(MMG5_pMesh mesh, MMG5_MIBOctree_s *q, MMG5_pMIBOctree *root, int64_t Z, int deepest) {

  if (q->is_filled == 0){
    MMG3D_split_MIBOctree_s ( mesh,q,&mesh->iboctree );
    q->is_filled = 1 ;
  }

  MMG5_MIBOctree_s* in_son = mesh->iboctree->root + q->sons[0];

  if (in_son->depth + 1 < deepest ){
    int64_t son_number = 0;
    int deplacement = 3*((*root)->depth_max - in_son->depth);
    son_number = (int64_t) 7 << deplacement;
    son_number = (Z & son_number) >> deplacement;

    MMG3D_split_in_balance(mesh,mesh->iboctree->root + q->sons[son_number], &mesh->iboctree, Z, deepest);
  }

  return 1;
}


/**
 * \param mesh pointer toward a mesh structure.
 * \param q pointer toward the MOctree cell
 * \param root pointer toward the octree root
 * \param deepest the maximal depth where unbalanced cells can exist
 *
 * \return 1 if success, 0 if fail.
 *
 * Recursive algorithme used to refine the octree for 2 : 1 balancing
 *
 */
static inline
int MMG3D_recursive_balance(MMG5_pMesh mesh, MMG5_MIBOctree_s *q, MMG5_pMIBOctree *root, int deepest) {

  int k = 0;

  if (q->depth == deepest ){
    int64_t Z_x;
    int64_t Z_y;
    int64_t Z_z;
    MMG5_MIBOctree_s* q_x;
    MMG5_MIBOctree_s* q_y;
    MMG5_MIBOctree_s* q_z;
    Z_x = MMG3D_external_neighbours(&mesh->iboctree, q->Z_coord, q->depth, 0);
    Z_y = MMG3D_external_neighbours(&mesh->iboctree, q->Z_coord, q->depth, 1);
    Z_z = MMG3D_external_neighbours(&mesh->iboctree, q->Z_coord, q->depth, 2);
    q_x = MMG3D_find_octree(mesh, &mesh->iboctree, Z_x);
    q_y = MMG3D_find_octree(mesh, &mesh->iboctree, Z_y);
    q_z = MMG3D_find_octree(mesh, &mesh->iboctree, Z_z);
    if (q_x->depth + 1 < deepest && Z_x >= 0)
      MMG3D_split_in_balance(mesh, q_x, &mesh->iboctree, Z_x, deepest);
    if (q_y->depth + 1 < deepest && Z_y >= 0)
      MMG3D_split_in_balance(mesh, q_y, &mesh->iboctree, Z_y, deepest);
    if (q_z->depth + 1 < deepest && Z_z >= 0)
      MMG3D_split_in_balance(mesh, q_z, &mesh->iboctree, Z_z, deepest);

  }else if (q->sons[0] != 0){
    for (k = 0 ; k < MMG3D_SIZE_OCTREESONS ; k++){
      MMG5_MIBOctree_s* son = mesh->iboctree->root + q->sons[k];
      MMG3D_recursive_balance(mesh, son, &mesh->iboctree, deepest);
    }
  }
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

  int depth = 0;
  int i = 0;
  for (depth = mesh->iboctree->max_depth_reached ; depth > 1 ; depth -= 1){
    for (i = 0 ; i < MMG3D_SIZE_OCTREESONS ; i ++)
      MMG3D_recursive_balance(mesh, mesh->iboctree->root + mesh->iboctree->root[0].sons[i], &mesh->iboctree, depth);
  }

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
  if ( !MMG5_scaleMesh_i( mesh, NULL,sol,MMG3D_IMMERSRAT ) ) {
    fprintf(stderr,"\n  ## Mesh scaling problem. Exit program.\n");
    return 0;
  }

  /* Surface mesh analysis (adjacency relation, tria orientation, normal at
   * points...) */
  if ( !MMGS_analys( mesh) ) {
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
  if ( !MMG3D_saveVTKMIBOctree(mesh,sol,"bbOctree.vtk") ) {
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
  // if ( !MMG3D_refine_octreeOnTria(mesh,sol) ) {
  //   fprintf(stderr,"\n  ## Octree refinement (over triangles) problem."
  //           " Exit program.\n");
  //   return 0;
  // }

  /* 2:1 balancing of the octree  */
  if ( !MMG3D_balance_octree(mesh,sol) ) {
    fprintf(stderr,"\n  ## Octree refinement (over triangles) problem."
            " Exit program.\n");
    return 0;
  }

  #ifndef NDEBUG
    // To remove when the octree creation will be ok
    if ( !MMG3D_saveVTKMIBOctree(mesh,sol,"finalOctree.vtk") ) {
      fprintf(stderr,"\n  ## Warning: unable to save the initial octree\n");
    }
  #endif

  /* Unscale the surface mesh and the octree */
  if ( !MMG5_unscaleMesh_i( mesh, NULL,sol,MMG3D_IMMERSRAT ) ) {
    fprintf(stderr,"\n  ## Mesh scaling problem. Exit program.\n");
    return 0;
  }
// #ifndef NDEBUG
//   // To remove when the octree creation will be ok
//   if ( !MMG3D_saveVTKMIBOctree(mesh,sol,"finalOctree.vtk") ) {
//     fprintf(stderr,"\n  ## Warning: unable to save the initial octree\n");
//   }
// #endif

  printf("vtk save was done\n");
  /* Memory release */
  MMG3D_delete_octree (mesh);

  return 1;
}
