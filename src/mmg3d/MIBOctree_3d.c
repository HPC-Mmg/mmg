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
 * \file MLSOctree_3d.c
 * \brief Tools to manage an octree mesh.
 * \author Gustavo Borges Valim (Enseirb-Matmeca)
 * \author Laurent Le (Enseirb-Matmeca)
 * \author Damien Sans (Enseirb-Matmeca)
 * \author Algiane Froehly (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Functions related to the creation of an octree and its subdivision on an
 * embedded grid
 *
 */

#include "mmg3d.h"

#define MMG3D_NOSONS 0

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MIBOctree
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and initialize a MIBOctree structure (octree for immersed boundary
 * capture).
 *
 */
int MMG3D_init_MIBOctree  ( MMG5_pMesh mesh, MMG5_pMIBOctree *q ) {
  size_t available;
  int    nbBitsInt,i;

  /** Root allocation */
  MMG5_ADD_MEM(mesh,sizeof(MMG5_MIBOctree),"MIBOctree root",
               return 0);
  MMG5_SAFE_MALLOC( *q,1, MMG5_MIBOctree, return 0);

  /** Computation of maximal depth allowing the computation of the z index */
  // nbBitsInt    = sizeof(int64_t)*8;
  // (*q)->depth_max = nbBitsInt - 1;
  // (*q)->depth_max = 20;      //value to use
  (*q)->depth_max = 6;    // test value


  /** Computation of initial maximal number of cells from the number of points */
  available  = mesh->memMax - mesh->memCur;
  // ncells_max may be underestimated: to be fitted
  available /=  sizeof(MMG5_MIBOctree_s);
  (*q)->ncells_max = 0;
  if ( MMG3D_SIZE_OCTREESONS*mesh->np < available ) {
    (*q)->ncells_max = MMG3D_SIZE_OCTREESONS*mesh->np;
    if ( (*q)->ncells_max +  MMG3D_SIZE_OCTREESONS*mesh->nt < available ) {
      (*q)->ncells_max +=  MMG3D_SIZE_OCTREESONS*mesh->nt;
    }
  }
  if ( !(*q)->ncells_max ) {
    (*q)->ncells_max = available;
  }
  if ( (*q)->ncells_max < MMG3D_SIZE_OCTREESONS ) {
    printf("  ## Error:%s:%d: Not enough memory to allocate octree array.\n",
           __func__,__LINE__);
    return 0;
  }

  // (*q)->ncells_max = 1000;
  /** Check that we have enough memory to allocate the array of cells */
  MMG5_ADD_MEM(mesh,(*q)->ncells_max*sizeof(MMG5_MIBOctree_s),"MIBOctree_s array",
               return 0);

  /** Array allocation */
  MMG5_SAFE_CALLOC( (*q)->root,(*q)->ncells_max, MMG5_MIBOctree_s, return 0);

  /** Creation of the 1 Cell of the minimal octree */
  (*q)->nxt   = 1;
  for ( i=0; i<MMG3D_SIZE_OCTREESONS; ++i ) {
    (*q)->root[0].sons[i] = MMG3D_NOSONS; //no sons
  }
  (*q)->root[0].depth   = 0;

  /** Creation of the linked list of cells */
  for ( i= 1; i<(*q)->ncells_max-1; ++i ) {
    (*q)->root[i].depth = i+1;
  }


  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param root pointer toward the octree root
 * \param depth depth of the cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Creation of a new MIBOctree cell at depth \a depth.
 *
 */
int MMG3D_newMIBOctree_s( MMG5_pMesh mesh,MMG5_pMIBOctree *root, int depth, int64_t z) {
  MMG5_MIBOctree_s* q;
  size_t            ncells_max;
  int               i;
  int               j;
  int               expnt;
  int64_t           a;


  if ( depth > (*root)->depth_max ) {
    printf(" ## Error:%s:%d: Attempt to create a new octree cell at depth %d"
           " deeper than depth_max (%d).\n",__func__,__LINE__,depth,(*root)->depth_max);
    return 0;
  }

  /* Reallocation of the cell array if too small */
  if ( !(*root)->nxt ) {
    printf(" ## Warning:%s:%d: Reallocation of the octree array."
           " We must pass here very few times.\n",__func__,__LINE__);
    ncells_max = mesh->gap*(*root)->ncells_max+1;
    MMG5_SAFE_REALLOC( (*root)->root,(*root)->ncells_max,ncells_max,MMG5_MIBOctree_s,"Octree array",
                       return 0);


    /* Linked list init */
    (*root)->nxt = (*root)->ncells_max;
    for ( i=(*root)->nxt; i<ncells_max-1; ++i ) {
      (*root)->root[i].depth = i+1;
    }
    (*root)->ncells_max = ncells_max;

  }

  i = (*root)->nxt;
  (*root)->nxt = (*root)->root[i].depth;
  q = (*root)->root + i;
  // printf("%d %d %d\n", i, (*root)->root[i].depth, q );
  j = i % MMG3D_SIZE_OCTREESONS;
  expnt = 3*((*root)->depth_max - depth);


  q->depth  = depth;

  q->is_filled = 0;
  q->Z_coord = z;

  a = 1;

  if ( j == 2){
    q->Z_coord += a << expnt;
  }else if (j == 3){
    q->Z_coord += a << (expnt + 1);
  }else if (j == 4){
    q->Z_coord += a << expnt;
    q->Z_coord += a << (expnt + 1);
  }else if (j == 5){
    q->Z_coord += a << (expnt + 2);
  }else if (j == 6){
    q->Z_coord += a << expnt;
    q->Z_coord += a << (expnt + 2);
  }else if (j == 7){
    q->Z_coord += a << (expnt + 1);
    q->Z_coord += a << (expnt + 2);
  }else if (j == 0){
    q->Z_coord += a << expnt;
    q->Z_coord += a << (expnt + 1);
    q->Z_coord += a << (expnt + 2);
  }
  // printf(" %d, son[%d], depth : %d , Z_coord : %lld \n", i, j, expnt/3, q->Z_coord );


  for ( i=0; i< MMG3D_SIZE_OCTREESONS; ++i ) {
    q->sons[i] = MMG3D_NOSONS; //no sons
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh structure
 * \param q pointer toward the MIBOctree cell
 * \param root pointer toward the octree root
 *
 * \return 1 if success, 0 if fail, -1 if we are at depth max.
 *
 * Split a MIBOctree cell \ref q into \ref  MMG3D_SIZE_OCTREESONS MLSOctree Cells.
 *
 */
int MMG3D_split_MIBOctree_s ( MMG5_pMesh mesh,MMG5_MIBOctree_s* q,MMG5_pMIBOctree *root) {
  int i;

  /* If the depth of the cell allow the subdivision */
  // printf("\n////// Father's Z_coord : %lld /////////////\n",q->Z_coord);
  if ( q->depth+1 < (*root)->depth_max) {
    for ( i=0; i <  MMG3D_SIZE_OCTREESONS; ++i ) {
      q->sons[i] = (*root)->nxt;

      if ( !MMG3D_newMIBOctree_s(mesh,root,q->depth + 1,q->Z_coord) ) {
        printf(" ## Error:%s:%d: Unable to split octree cell\n",
               __func__,__LINE__);
        return 0;
      }
    }
    // q->is_filled = 0;

  }
  else {
    return -1;
  }

  return 1;
}

/**
 * \param mesh toward the mesh structure
 * \param q pointer toward the MIBOctree
 *
 * \return 1 if success, 0 if fail.
 *
 * Free a MIBOctree root.
 *
 */
int MMG3D_free_MIBOctree  ( MMG5_pMesh mesh,MMG5_pMIBOctree* q ) {

  MMG5_DEL_MEM(mesh,(*q)->root);

  MMG5_DEL_MEM(mesh,(*q));

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param id_cell index of the cell to delete
 * \param root pointer toward the octree root
 *
 * \return 1 if success, 0 if fail.
 *
 * Suppression of a new MIBOctree cell.
 *
 */
int MMG3D_delMIBOctree_s( MMG5_pMesh mesh,int id_cell,MMG5_pMIBOctree *root ) {
  MMG5_MIBOctree_s* q;

  q = (*root)->root + id_cell;

  q->depth     = (*root)->nxt;
  (*root)->nxt = id_cell;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param q pointer toward the MIBOctree cell
 * \param root pointer toward the octree root
 *
 * \return 1 if success, 0 if fail, -1 if we are at depth max.
 *
 * Merge the sons of an octree cell
 *
 */
int  MMG3D_merge_MIBOctree_s ( MMG5_pMesh mesh,MMG5_MIBOctree_s* q,MMG5_pMIBOctree *root) {
  int i,is_filled;

  is_filled = 0;

  if ( q->sons[0] == MMG3D_NOSONS ) {
    printf(" ## Error:%s:%d: Unable to merge a leaf cell\n",
           __func__,__LINE__);
  }

  for ( i=0; i <  MMG3D_SIZE_OCTREESONS; ++i ) {
    if ( !is_filled && (*root)->root[q->sons[i]].is_filled ) {
      is_filled = 1;
    }

    if ( !MMG3D_delMIBOctree_s(mesh,q->sons[i],root) ) {
      printf(" ## Error:%s:%d: Unable to merge octree cell\n",
             __func__,__LINE__);
      return 0;
    }
    q->sons[i] = MMG3D_NOSONS;
  }
  q->is_filled = is_filled;

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure
 * \param q pointer toward the MIBOctree cell
 * \param root pointer toward the octree root
 * \param np number of points (to be updated)
 * \param nc number of leaf cells (to be updated)
 *
 * Count the number of cells and the number of corners of leaf cells in the
 * octree branch.
 *
 */
void MMG3D_count_MIBOctreeEntities(MMG5_MIBOctree_s* q,MMG5_pMIBOctree root,int *np,int *nc) {
  int k;


  if ( q->sons[0] == MMG3D_NOSONS ) {
    (*np) += 8;
    ++(*nc);

  }
  else {
    for ( k=0; k< MMG3D_SIZE_OCTREESONS; ++k ) {
      MMG3D_count_MIBOctreeEntities(root->root + q->sons[k],root,np,nc);
    }
  }

  return;
}

/**
 * \param inm file pointer
 * \param mesh pointer toward the mesh structure
 * \param q pointer toward the MIBOctree cell
 * \param root pointer toward the octree root
 * \param x_min minimal coordinates of the branch in the x direction
 * \param y_min minimal coordinates of the branch in the y direction
 * \param z_min minimal coordinates of the branch in the z direction
 * \param x_max maximal coordinates of the branch in the x direction
 * \param y_max maximal coordinates of the branch in the y direction
 * \param z_max maximal coordinates of the branch in the z direction
 *
 * Write the Corners of leaf cells in the octree branch.
 *
 */
void MMG3D_write_MIBOctreeCoor(FILE *inm,MMG5_MIBOctree_s* q,MMG5_pMIBOctree root,
                               double x_min,double y_min,double z_min,
                               double x_max,double y_max,double z_max) {
  double x_mid,y_mid,z_mid;

  if ( q->sons[0] == MMG3D_NOSONS ) {
    /* Bottom left front corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_min,y_min,z_min);
    /* Bottom right front corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_max,y_min,z_min);
    /* Bottom right back corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_max,y_max,z_min);
    /* Bottom left back corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_min,y_max,z_min);
    /* Top left front corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_min,y_min,z_max);
    /* Top right front corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_max,y_min,z_max);
    /* Top right back corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_max,y_max,z_max);
    /* Top left back corner */
    fprintf(inm,"%.15lg %.15lg %.15lg\n",x_min,y_max,z_max);
  }
  else {
    x_mid = ( x_max + x_min ) / 2.;
    y_mid = ( y_max + y_min ) / 2.;
    z_mid = ( z_max + z_min ) / 2.;

    /* Bottom left front branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[0],root,
                              x_min,y_min,z_min,x_mid,y_mid,z_mid);

    /* Bottom right front branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[1],root,
                              x_mid,y_min,z_min,x_max,y_mid,z_mid);

    /* Bottom left back branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[2],root,
                              x_min,y_mid,z_min,x_mid,y_max,z_mid);

    /* Bottom right back branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[3],root,
                              x_mid,y_mid,z_min,x_max,y_max,z_mid);

    /* Top left front branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[4],root,
                              x_min,y_min,z_mid,x_mid,y_mid,z_max);

    /* Top right front branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[5],root,
                              x_mid,y_min,z_mid,x_max,y_mid,z_max);

    /* Top left back branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[6],root,
                              x_min,y_mid,z_mid,x_mid,y_max,z_max);

    /* Top right back branch */
    MMG3D_write_MIBOctreeCoor(inm,root->root + q->sons[7],root,
                              x_mid,y_mid,z_mid,x_max,y_max,z_max);



  }

  return;
}

/**
 * \param inm file pointer
 * \param mesh pointer toward the mesh structure
 * \param q pointer toward the MIBOctree cell
 * \param root pointer toward the octree root
 * \param np index of current point (to be updated)
 *
 * Write the connectivity of the leaves cells in the octree branch.
 *
 */
void MMG3D_write_MIBOctreeVertices(FILE *inm,MMG5_MIBOctree_s* q,MMG5_pMIBOctree root,int *np) {
  int k;


  if ( q->sons[0] == MMG3D_NOSONS ) {
    fprintf(inm,"%d %d %d %d %d %d %d %d %d\n",MMG3D_SIZE_OCTREESONS,
            (*np),(*np)+1,(*np)+2,(*np)+3,(*np)+4,(*np)+5,(*np)+6,(*np)+7);
    (*np) += 8;
  }
  else {
    for ( k=0; k< MMG3D_SIZE_OCTREESONS; ++k ) {
      MMG3D_write_MIBOctreeVertices(inm,root->root + q->sons[k],root,np);
    }
  }

  return;
}
