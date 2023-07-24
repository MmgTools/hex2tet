/**
 * \file libhex2tet.c
 * \brief main library function: cut an hexahedral mesh into a tetrahedral one
 *
 * \author Cecile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

int H2T_libhex2tet(MMG5_pMesh mmgMesh,int** hexa,MMG5_int nbhexa) {
  Hedge             hed2;
  int               norient;
  MMG5_int          *adjahex,k,ier;
  
  if ( mmgMesh->info.imprim ) {
    fprintf(stdout,"\n  -- HEX2TET, Release %s (%s) \n",H2T_VER,H2T_REL);
    fprintf(stdout,"     %s\n",H2T_CPY);
    fprintf(stdout,"     %s %s\n",__DATE__,__TIME__);
  }

  /* chk orientation */
  norient = H2T_chkorient(mmgMesh,*hexa,nbhexa);
  if ( norient )
    fprintf(stdout,"\n  -- WARNING: %8d HEXA REORIENTED\n",norient);

  /* hexa adjacency */
  adjahex = NULL;
  adjahex = (MMG5_int*)calloc(6*nbhexa+7,sizeof(long long));
  assert(adjahex);

  if(!H2T_hashHexa(*hexa,adjahex,nbhexa)) return H2T_STRONGFAILURE;

  /* cut hexa into tet */
  hed2.size  = 6*nbhexa;
  hed2.hnxt  = 6*nbhexa;
  //hed2.nhmax = (MMG5_int)(16*6*nbhexa);
  hed2.nhmax = (MMG5_int)(2*6*nbhexa);
  hed2.item  = NULL;
  hed2.item  = (hedge*)calloc(hed2.nhmax+1,sizeof(hedge));

  for (k=6*nbhexa; k<hed2.nhmax; k++) {
    hed2.item[k].nxt = k+1;
  }

  if ( 6*nbhexa >= mmgMesh->nemax ) {
    printf("\n  -- ERROR: Not enough memory to store the final mesh"
           " (max number of tetra=%" MMG5_PRId " while the number of hexa times"
           " 6 is % " MMG5_PRId ".\n",
           mmgMesh->nemax,nbhexa*6);
    return H2T_STRONGFAILURE;
  }

  ier = H2T_cuthex(mmgMesh,&hed2,*hexa,adjahex,nbhexa);

  if ( !ier )
    ier = H2T_STRONGFAILURE;
  else
    ier = H2T_SUCCESS;

  free(adjahex);
  free(hed2.item);

  return ier;
}
