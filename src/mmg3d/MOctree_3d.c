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
 * \file MOctree_3d.c
 * \brief Tools to manage an octree mesh.
 * \author Juliette Busquet (Enseirb-Matmeca)
 * \author Antoine Huc (Enseirb-Matmeca)
 * \author Fanny Kuhn (Enseirb-Matmeca)
 * \author Algiane Froehly (Inria)
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 * Functions related to the creation of an octree from a grid and the associated
 * signed distance function. The octree cells are allowed to merge if they do
 * not cross the 0 levels-set of the signed distance function. At the end, the
 * octree must be 2:1 balanced (only 1 depth level of difference betwwen to
 * adjacent cells).
 *
 */

#include "mmg3d.h"

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree root
 * \param ip index of the bottom-left-front corner of the root cell
 * \param length length of the octree in each direction
 * \param depth_max maximal depth of the octree
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree structure.
 *
 */
int MMG3D_init_MOctree  ( MMG5_pMesh mesh, MMG5_pMOctree *q, int ip,
                          double length[3],int depth_max ) {


  /** Root allocation */
  MMG5_ADD_MEM(mesh,sizeof(MMG5_MOctree),"MOctree root",
               return 0);
  MMG5_SAFE_MALLOC( *q,1, MMG5_MOctree, return 0);

  (*q)->length[0] = length[0];
  (*q)->length[1] = length[1];
  (*q)->length[2] = length[2];

  (*q)->depth_max = depth_max;

  (*q)->nspan_at_root = pow(2,depth_max);

  /** Check that we have enough memory to allocate a new cell */
  MMG5_ADD_MEM(mesh,sizeof(MMG5_MOctree_s),"initial MOctree cell",
               return 0);

  /** New cell allocation */
  MMG5_SAFE_MALLOC( (*q)->root,1, MMG5_MOctree_s, return 0);
  MMG3D_init_MOctree_s( mesh,(*q)->root,ip,0,0);
  (*q)->root->nsons = 0;

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 * \param ip index of the bottom-left-front corner of the cell
 * \param depth cell's depth
 * \param split_ls 1 if the cell is intersected by the level-set
 *
 * \return 1 if success, 0 if fail.
 *
 * Allocate and init an MOctree Cell
 *
 */
int MMG3D_init_MOctree_s( MMG5_pMesh mesh, MMG5_MOctree_s* q,int ip, int depth,int8_t split_ls ) {

  q->father = NULL;
  q->sons   = NULL;
  q->nsons  = 0;

  q->blf_ip = ip;
  q->depth  = depth;

  q->split_ls = split_ls;
  q->leaf = 0;
  q->coordoct[0]=0;
  q->coordoct[1]=0;
  q->coordoct[2]=0;
  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Set the parameter split_ls of a cell to 1 if the level-set intersect the cell for every cell of the octree max.
 *
 */
int  MMG3D_set_splitls_MOctree ( MMG5_pMesh mesh, MMG5_MOctree_s* q, MMG5_pSol sol) {
  int ip;
  int FDL,FDR,BDL,BDR,FUL,FUR,BUL,BUR;// ip of the 8 vertices of an octree cell
  ip=q->blf_ip;

  FDL=ip; // Front Down Left corner
  FDR=FDL+1;
  BDL=ip+mesh->freeint[0];
  BDR=BDL+1;
  FUL=ip+mesh->freeint[0]*mesh->freeint[1];
  FUR=FUL+1;
  BUL=ip+mesh->freeint[0]*(1+mesh->freeint[1]);
  BUR=BUL+1;

  if(sol->m[FDL]*sol->m[FDR]<0 || sol->m[FDR]*sol->m[BDL]<0 ||\
     sol->m[BDL]*sol->m[BDR]<0 || sol->m[BDR]*sol->m[FUL]<0 ||\
     sol->m[FUL]*sol->m[FUR]<0 || sol->m[FUR]*sol->m[BUL]<0 || (sol->m[BUL]*sol->m[BUR]<0))
  {
    q->split_ls=1;
  }
  else
  {
    q->split_ls=0;
  }

  return 1;
}

/**
 * \param mesh pointer toward a MMG5 mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Split an MOctree cell \ref q into \ref MMG3D_SIZE_OCTREESONS MOctree Cells.
 * \ref q must be a leaf.
 *
 */
 int  MMG3D_split_MOctree_s ( MMG5_pMesh mesh,MMG5_MOctree_s* q,MMG5_pSol sol, int* span) {
   int ip,i,power,blf_ip_first_son;

   int depth_max = mesh->octree->depth_max;

   if(q->depth < depth_max)
   {
     MMG5_SAFE_MALLOC(q->sons,q->nsons, MMG5_MOctree_s, return 0);
     for(i=0; i<q->nsons; i++)
     {
       MMG3D_init_MOctree_s(mesh, &q->sons[i], 0, q->depth + 1, 0);
       //calculus of octree coordinates and ip
       int span_y,span_z;
       int power = pow(2,depth_max-(q->depth+1));
       //int ncells_x  = (mesh->freeint[0]-1)/power + 1;
       int ncells_x = pow(2,q->depth+1) + 1;
       int ncells_y = pow(2,q->depth+1) + 1;
       int ncells_z = pow(2,q->depth+1) + 1;
       if(ncells_x > mesh->freeint[0])
       {
         ncells_x = mesh->freeint[0];
       }
       if(ncells_y > mesh->freeint[1])
       {
         ncells_y = mesh->freeint[1];
       }
       if(ncells_z > mesh->freeint[2])
       {
         ncells_z = mesh->freeint[2];
       }

       int ncells_xy = ncells_x * ncells_y;

       span_y = ncells_x;
       span_z = ncells_xy;

       q->sons[i].coordoct[0]=q->coordoct[0];
       q->sons[i].coordoct[1]=q->coordoct[1];
       q->sons[i].coordoct[2]=q->coordoct[2];

       if(i==0)
       {
         q->sons[i].blf_ip = q->coordoct[0]/power + q->coordoct[1]/power*span_y + q->coordoct[2]/power*span_z + 1;
       }
       else if(i==1)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + 1;
         q->sons[i].coordoct[0]+=power;
       }
       else if(i==2)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_y;
         q->sons[i].coordoct[1]+=power;
       }
       else if(i==3)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_y + 1;
         q->sons[i].coordoct[0]+=power;
         q->sons[i].coordoct[1]+=power;
       }
       else if(i==4)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_z;
         q->sons[i].coordoct[2]+=power;
       }
       else if(i==5)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_z + 1;
         q->sons[i].coordoct[0]+=power;
         q->sons[i].coordoct[2]+=power;
       }
       else if(i==6)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_z + span_y;
         q->sons[i].coordoct[1]+=power;
         q->sons[i].coordoct[2]+=power;
       }
       else if(i==7)
       {
         q->sons[i].blf_ip = q->sons[0].blf_ip + span_z + span_y + 1;
         q->sons[i].coordoct[0]+=power;
         q->sons[i].coordoct[1]+=power;
         q->sons[i].coordoct[2]+=power;
       }

       q->sons[i].father = q;
       q->sons[i].nsons = 8;
       MMG3D_split_MOctree_s(mesh, &q->sons[i], sol, span);
     }
   }
   else{
     /*calculus of ip for leaves*/
     //q->blf_ip=q->coordoct[2]*(pow(2,2*depth_max)+1)+q->coordoct[1]*(pow(2,depth_max)+1)+q->coordoct[0]+1;
     //moins couteux ?
     //q->sons[i].blf_ip=q->sons[i].coordoct[2]*dx*dy+q->sons[i].coordoct[1]*dx+q->sons[i].coordoct[0]+1;
     q->leaf=1;
     MMG3D_set_splitls_MOctree (mesh, q, sol);
   }

   return 1;
 }

/**
 * \param q pointer toward the MOctree
 *
 * \return 1 if success, 0 if fail.
 *
 * Free a MOctree.
 *
 */
int MMG3D_free_MOctree  ( MMG5_pMOctree** q, MMG5_pMesh mesh) {
  MMG5_DEL_MEM(mesh,*q);
  return 1;
}

/**
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Free a MOctree cell.
 *
 */
int MMG3D_free_MOctree_s( MMG5_MOctree_s* q, MMG5_pMesh mesh) {
  MMG5_DEL_MEM(mesh,q);
  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param q pointer toward the MOctree cell
 *
 * \return 1 if success, 0 if fail.
 *
 * Merge the \ref q MOctree cell (remove the \ref MMG3D_SIZE OCTREESONS sons of
 * \ref q, all sons being leafs).
 *
 */
int  MMG3D_merge_MOctree_s ( MMG5_MOctree_s* q, MMG5_pMesh mesh) {

  q->nsons = 0;
  q->leaf = 1;
  q->blf_ip=q->sons[0].blf_ip;//récupérer la ls à l'origine de la cellule
  MMG3D_free_MOctree_s(q->sons, mesh);

  return 1;
}

/**
 * \param mesh pointer toward the mesh
 * \param q pointer toward the MOctree cell
 * \param span span between a corner and the other corner in one given direction
 * \param ip0 pointer toward the index of the bottom left front corner of the cell (to fill)
 * \param ip1 pointer toward the index of the bottom right front corner of the cell (to fill)
 * \param ip2 pointer toward the index of the bottom right back corner of the cell (to fill)
 * \param ip3 pointer toward the index of the bottom left back corner of the cell (to fill)
 * \param ip4 pointer toward the index of the top left front corner of the cell (to fill)
 * \param ip5 pointer toward the index of the top right front corner of the cell (to fill)
 * \param ip6 pointer toward the index of the top right back corner of the cell (to fill)
 * \param ip7 pointer toward the index of the top left back corner of the cell (to fill)
 *
 * \return 1 if success, 0 if fail.
 *
 * Compute the indices of the corners of the octree cell (vertices are numbering
 * as in the initial grid)
 *
 */
int MMG3D_get_MOctreeCornerIndices ( MMG5_pMesh mesh, MMG5_MOctree_s *q,int span,int *ip0,
                                     int *ip1, int *ip2,int *ip3,int *ip4,int *ip5,
                                     int *ip6,int *ip7 ) {

  int span_y,span_z;
  int ncells_x  = mesh->freeint[0];
  int ncells_xy = mesh->freeint[0]*mesh->freeint[1];

  span_y = ncells_x;
  span_z = ncells_xy;

  *ip0 = q->blf_ip;
  *ip1 = q->blf_ip + 1 ;
  *ip2 = q->blf_ip + span_y + 1;
  *ip3 = q->blf_ip + span_y;
  *ip4 = q->blf_ip + span_z;
  *ip5 = q->blf_ip + span_z + 1;
  *ip6 = q->blf_ip + span_z + span_y + 1;
  *ip7 = q->blf_ip + span_z + span_y;

  return 1;
}
/**
 * \param mesh pointer toward the mesh
 * \param q pointer toward the MOctree cell
 * \param span span between a corner and the other corner in one given direction
 * \param np pointer toward the number of used points (to fill)
 * \param nc pointer toward the number of cell leafs (to fill)
 *
 * \return 1 if success, 0 if fail.
 *
 * Mark as used the points that are at the corners of the octree cells
 * leafs. Count the number of use points and the number of leafs.
 *
 */
int  MMG3D_mark_MOctreeCellCorners ( MMG5_pMesh mesh, MMG5_MOctree_s* q,int *span,int *np,int *nc ) {
  MMG5_pPoint ppt;
  int         i,ip[8];

  if ( q->leaf ) {

    if ( !MMG3D_get_MOctreeCornerIndices ( mesh,q,*span,ip,ip+1,ip+2,ip+3,ip+4,ip+5,ip+6,ip+7 ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute the indices of the"
              " corners of the octree cell.\n",__func__);
      return 0;
    }

    for ( i=0; i<8; ++i ) {
      ppt = &mesh->point[ip[i]];
      if ( !MG_VOK(ppt) ) {
        ++(*np);
        ppt->tag &= ~MG_NUL;
      }
    }
    ++(*nc);
  }
  else {
    (*span) /= 2;

    for ( i=0; i<q->nsons; ++i ) {
      if ( !MMG3D_mark_MOctreeCellCorners ( mesh,&q->sons[i],span,np,nc ) ) return 0;
    }
  }

  return 1;
}


/**
 * \param mesh pointer toward the mesh
 * \param q pointer toward the MOctree cell
 * \param span span between a corner and the other corner in one given direction
 * \param inm pointer toward the file in which we save the octree cells
 *
 * \return 1 if success, 0 if fail.
 *
 * Write the hexahedron associated to the octree leafs.
 *
 */
int  MMG3D_write_MOctreeCell ( MMG5_pMesh mesh, MMG5_MOctree_s* q,int *span,FILE *inm ) {
  int i,ip0,ip1,ip2,ip3,ip4,ip5,ip6,ip7;
  static int nvert = 8;

  if ( q->leaf ) {

    if ( !MMG3D_get_MOctreeCornerIndices ( mesh,q,*span,&ip0,&ip1,&ip2,&ip3,
                                           &ip4,&ip5,&ip6,&ip7 ) ) {
      fprintf(stderr,"\n  ## Error: %s: unable to compute the indices of the"
              " corners of the octree cell.\n",__func__);
      return 0;
    }

    fprintf(inm,"%d %d %d %d %d %d %d %d %d\n",nvert,mesh->point[ip0].tmp,
            mesh->point[ip1].tmp,mesh->point[ip2].tmp,mesh->point[ip3].tmp,
            mesh->point[ip4].tmp,mesh->point[ip5].tmp,mesh->point[ip6].tmp,
            mesh->point[ip7].tmp);
  }
  else {
    (*span) /= 2;

    for ( i=0; i<q->nsons; ++i ) {
      if ( !MMG3D_write_MOctreeCell ( mesh,&q->sons[i],span,inm ) ) {
        return 0;
      }
    }
  }

  return 1;
}

int MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(MMG5_pMesh mesh, MMG5_MOctree_s* q, char Direction[20], MMG5_MOctree_s* Neighbour)
{
  int i;
  MMG5_MOctree_s* Temp_Neighbour;

  if(Direction == "UP")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==1 || i==2 || i==3)  && q==&q->father->sons[i]) // Is q a south child child?
      {
        Neighbour=&q->father->sons[i+4];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }

    //q is guaranted to be a North child
    for (i=0; i<8; i++)
    {
      if((i==4 || i==5 || i==6 || i==7) && q==&q->father->sons[i]) // Is q a south child child?
      {
        Neighbour=&Temp_Neighbour->father->sons[i-4];
        return 1;
      }
    }
  }

  if(Direction == "DOWN")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==4 || i==5 || i==6 || i==7)  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-4];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }


    for (i=0; i<8; i++)
    {
      if((i==0 || i==1 || i==2 || i==3) && q==&q->father->sons[i])
      {
        Neighbour=&Temp_Neighbour->father->sons[i+4];
        return 1;
      }
    }
  }

  if(Direction == "LEFT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==1 || i==3 || i==5 || i==7)  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-1];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }


    for (i=0; i<8; i++)
    {
      if((i==0 || i==2 || i==4 || i==6) && q==&q->father->sons[i])
      {
        Neighbour=&Temp_Neighbour->father->sons[i+1];
        return 1;
      }
    }
  }


  if(Direction == "RIGHT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==2 || i==4 || i==6) && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+1];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }


    for (i=0; i<8; i++)
    {
      if((i==1 || i==3 || i==5 || i==7)  && q==&q->father->sons[i])
      {
        Neighbour=&Temp_Neighbour->father->sons[i-1];
        return 1;
      }
    }
  }

  if(Direction == "BACK")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==1 || i==4 || i==5) && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+2];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }


    for (i=0; i<8; i++)
    {
      if((i==2 || i==3 || i==6 || i==7)  && q==&q->father->sons[i])
      {
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
    }
  }



  if(Direction == "FRONT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==2 || i==3 || i==6 || i==7)  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-2];
        return 1;
      }
    }

    MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,Direction,Temp_Neighbour);
    if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
    {
      Neighbour=Temp_Neighbour;
      return 1;
    }


    for (i=0; i<8; i++)
    {
      if((i==0 || i==1 || i==4 || i==5) && q==&q->father->sons[i])
      {
        Neighbour=&Temp_Neighbour->father->sons[i+2];
        return 1;
      }
    }
  }



  if(Direction == "UPLEFT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==1 || i==3 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+3];
      }
      else if((i==0 || i==2 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+5];
        return 1;
      }
      else if((i==5 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UP",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-5];
        return 1;
      }
      else if((i==4 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPLEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
    }
  }

  if(Direction == "UPRIGHT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==2 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+5];
      }
      else if((i==1 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+3];
        return 1;
      }
      else if((i==4 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UP",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-3];
        return 1;
      }
      else if((i==5 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPRIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-5];
        return 1;
      }
    }
  }

  if(Direction == "DOWNLEFT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==5 || i==7 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-5];
      }
      else if((i==4 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-3];
        return 1;
      }
      else if((i==1 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UP",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+3];
        return 1;
      }
      else if((i==0 || i==2 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"DOWNLEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-5];
        return 1;
      }
    }
  }


  if(Direction == "DOWNRIGHT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==4 || i==6 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-3];
      }
      else if((i==5 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-5];
        return 1;
      }
      else if((i==0 || i==2 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"DOWN",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+5];
        return 1;
      }
      else if((i==1 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"DOWNRIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
    }
  }


  if(Direction == "UPFRONT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==2 || i==3 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+2];
      }
      else if((i==0 || i==1 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"FRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+6];
        return 1;
      }
      else if((i==6 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UP",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-6];
        return 1;
      }
      else if((i==4 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPFRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
    }
  }


  if(Direction == "UPBACK")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==1 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+6];
      }
      else if((i==2 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"BACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+2];
        return 1;
      }
      else if((i==4 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UP",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
      else if((i==6 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPBACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-6];
        return 1;
      }
    }
  }



  if(Direction == "DOWNFRONT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==6 || i==7 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-6];
      }
      else if((i==4 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"FRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-2];
        return 1;
      }
      else if((i==2 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"DOWN",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+2];
        return 1;
      }
      else if((i==0 || i==1 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPFRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+6];
        return 1;
      }
    }
  }


  if(Direction == "DOWNBACK")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==4 || i==5 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-2];
      }
      else if((i==6 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"BACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-6];
        return 1;
      }
      else if((i==0 || i==1 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"DOWN",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+6];
        return 1;
      }
      else if((i==2 || i==3 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"UPFRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+2];
        return 1;
      }
    }
  }


  if(Direction == "LEFTFRONT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==3 || i==7 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-3];
      }
      else if((i==1 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"FRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+1];
        return 1;
      }
      else if((i==2 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-1];
        return 1;
      }
      else if((i==0 || i==4 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFTFRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+3];
        return 1;
      }
    }
  }


  if(Direction == "LEFTBACK")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==1 || i==5 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+1];
      }
      else if((i==3 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"BACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-3];
        return 1;
      }
      else if((i==0 || i==4 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+3];
        return 1;
      }
      else if((i==2 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"LEFTBACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-1];
        return 1;
      }
    }
  }


  if(Direction == "RIGHTFRONT")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==2 || i==6 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i-1];
      }
      else if((i==0 || i==4 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"FRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+3];
        return 1;
      }
      else if((i==3 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-3];
        return 1;
      }
      else if((i==1 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHTFRONT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+1];
        return 1;
      }
    }
  }


  if(Direction == "RIGHTBACK")
  {
    if (q==mesh->octree->root) // q is root = q has no neighbour
    {
      Neighbour=NULL;
      return 1;
    }

    for (i=0; i<8; i++)
    {
      if((i==0 || i==4 )  && q==&q->father->sons[i])
      {
        Neighbour=&q->father->sons[i+3];
      }
      else if((i==2 || i==6 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"BACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-1];
        return 1;
      }
      else if((i==1 || i==5 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHT",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i+1];
        return 1;
      }
      else if((i==3 || i==7 )  && q==&q->father->sons[i])
      {
        MMG3D_find_Neighbour_of_Bigger_or_Equal_Size(mesh,q->father,"RIGHTBACK",Temp_Neighbour);
        if(Temp_Neighbour==NULL || Temp_Neighbour->leaf==1)
        {
          Neighbour=Temp_Neighbour;
          return 1;
        }
        Neighbour=&Temp_Neighbour->father->sons[i-3];
        return 1;
      }
    }
  }
}
