/**
 * Example of call from a C code of the hex2tet library (convert an
 * array of hexa into a MMG5 tetrahedral mesh)
 *
 * \author Vivien Pianet
 * \author Algiane Froehly (InriaSoft)
 * \version 1
 * \copyright GNU Lesser General Public License.
 */


#include "hex2tet/libhex2tet.h"

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <string.h>


int main ( int argc, char* argv[] ) {

  const int Width  = 50;
  const int Height = 20;
  const int Depth  = 10;

  MMG5_pMesh mmgMesh = NULL;
  MMG5_pSol  mmgSol  = NULL;

  int        ier;
  char       *fileout;

  fprintf(stdout,"  -- HEX2TET EXAMPLE \n");

  if ( argc != 2 ) {
    printf(" Usage: %s fileout \n",argv[0]);
    return 1;
  }

  /* Name and path of the output mesh file */
  fileout = (char *) calloc(strlen(argv[1]) + 1, sizeof(char));
  if ( fileout == NULL ) {
    perror("  ## Memory problem: calloc");
    exit(EXIT_FAILURE);
  }
  strcpy(fileout,argv[1]);

  /** Step 1: MMG5 mesh allocation */
  int nbVertices = Width*Height*Depth;
  int nbHex = (Width - 1)    * (Height - 1) * (Depth-1);
  int nbQuad = 2*(Width-1)*(Height-1)+2*(Height-1)*(Depth-1)+2*(Width-1)*(Depth-1);

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /** Step 2: Set the mesh size by giving the number of vertices and
   * the number of hexa of the hexahedral mesh */
  if ( H2T_Set_meshSize(mmgMesh, nbVertices, nbHex, nbQuad, 0) != 1)
    exit(EXIT_FAILURE);

  /** Step 3: Give the mesh vertices to Hex2tet (hexahedra vertices) */
  int imax = Width;
  int jmax = Height;
  int kmax = Depth;

  int ref = 0;

  for ( int k=0; k<kmax; ++k ) {
    for ( int j=0; j<jmax; ++j ) {
      for ( int i=0; i<imax; ++i ) {
        /* Vertices insertion */
        int pointNumber = (k*jmax*imax) + (j*imax) + i + 1;
        double v[3];
        v[0] = (double)i;
        v[1] = (double)j;
        v[2] = (double)k;
        if ( H2T_Set_vertex(mmgMesh, v[0], v[1], v[2], ref, pointNumber) != 1)
          exit ( EXIT_FAILURE );
      }
    }
  }

  /** Step 4: Fill the array of the hexa connectivity (hexTab) such as:
   * hexTab [ 9*k + i] is the vertex number i of the k^th hexahedra
   * (for k from 1 to nbHex and i from 0 to 7) and  hexTab [ 9*k + 8] is the
   * reference of the k^th hexa */
  int *hexTab;
  hexTab = (int*)malloc(9*(nbHex+1)*sizeof(int));

  for ( int k=0; k<kmax-1; ++k ) {
    for ( int j=0; j<jmax-1; ++j ) {
      for ( int i=0; i<imax-1; ++i ) {
        int hexaNumber     = (k*(jmax-1)*(imax-1)) + (j*(imax-1)) + i + 1;
        /* Hexahedra definition and storage
         * f = front b = back || u = up vs. d = down || l = left vs. r = right */

        int fdl,fdr,bdl,bdr,ful,fur,bul,bur;
        fdl = (k * imax * jmax) + (j * imax) + i + 1;
        fdr = fdl + 1;
        bdl = ((k+1) * imax * jmax) + (j * imax) + i + 1;
        bdr = bdl + 1;
        ful = (k * imax * jmax) + ((j+1) * imax) + i + 1;
        fur = ful + 1;
        bul = ((k+1) * imax * jmax) + ((j+1) * imax) + i + 1;
        bur = bul + 1;

        if( H2T_Set_hexahedron(hexTab,fdl,fdr,bdr,bdl,ful,fur,bur,bul,
                               ref,hexaNumber) != 1 )
          exit ( EXIT_FAILURE );
      }
    }
  }


  /** Set boundary quadrangles */
  int quadNumber = 0;

  /* Lower face */
  ref = 1;
  int j = 0;
  for ( int k=0; k<kmax-1; ++k ) {
    for ( int i=0; i<imax-1; ++i ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = ((k+1) * imax * jmax) + (j * imax) + i + 1;
      dr = dl + 1;
      ul = (k * imax * jmax) + (j * imax) + i + 1;
      ur = ul + 1;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }

  /* Upper face */
  ref = 5;
  j = jmax-1;
  for ( int k=0; k<kmax-1; ++k ) {
    for ( int i=0; i<imax-1; ++i ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = ((k+1) * imax * jmax) + (j * imax) + i + 1;
      dr = dl + 1;
      ul = (k * imax * jmax) + (j * imax) + i + 1;
      ur = ul + 1;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }

  /* Front face */
  ref = 2;
  int k = 0;
  for ( int j=0; j<jmax-1; ++j ) {
    for ( int i=0; i<imax-1; ++i ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = (k * imax * jmax) + (j * imax) + i + 1;
      dr = dl + 1;
      ul = (k * imax * jmax) + ((j+1) * imax) + i + 1;
      ur = ul + 1;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }

  /* Back face */
  ref = 4;
  k = kmax-1;
  for ( int j=0; j<jmax-1; ++j ) {
    for ( int i=0; i<imax-1; ++i ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = (k * imax * jmax) + (j * imax) + i + 1;
      dr = dl + 1;
      ul = (k * imax * jmax) + ((j+1) * imax) + i + 1;
      ur = ul + 1;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }

  /* Right face */
  ref = 3;
  int i = imax-1;
  for ( int k=0; k<kmax-1; ++k ) {
    for ( int j=0; j<jmax-1; ++j ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = (k * imax * jmax) + (j * imax) + i + 1;
      dr = dl + imax*jmax;
      ul = (k * imax * jmax) + ((j+1) * imax) + i + 1;
      ur = ul + imax*jmax;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }

  /* Left face */
  ref = 5;
  i = 0;
  for ( int k=0; k<kmax-1; ++k ) {
    for ( int j=0; j<jmax-1; ++j ) {
      ++quadNumber;
      /* Quadrilateral definition and storage
       * u = up vs. d = down || l = left vs. r = right */

      int dl,dr,ul,ur;
      dl = (k * imax * jmax) + (j * imax) + i + 1;
      dr = dl + imax*jmax;
      ul = (k * imax * jmax) + ((j+1) * imax) + i + 1;
      ur = ul + imax*jmax;

      if( H2T_Set_quadrilateral(mmgMesh,dl,dr,ur,ul,ref,quadNumber) != 1 )
        exit ( EXIT_FAILURE );
    }
  }
  assert(quadNumber == mmgMesh->nquad);


//  /* Save he hexa mesh */
//  FILE* inm = fopen("hexa.mesh","w");
//
//  fprintf(inm,"MeshVersionFormatted 2\n");
//  fprintf(inm,"\n\nDimension 3\n");
//  fprintf(inm,"\n\nVertices\n");
//
//  fprintf(inm,"%d\n",mmgMesh->np);
//
//  for (int k=1; k<=mmgMesh->np; k++) {
//    MMG5_pPoint ppt = &mmgMesh->point[k];
//    fprintf(inm,"%.15lg %.15lg %.15lg %d\n",ppt->c[0],ppt->c[1],ppt->c[2],abs(ppt->ref));
//  }
//
//  fprintf(inm,"\n\nHexahedra %d\n",nbHex);
//  for (int k=1; k<=nbHex; k++) {
//    fprintf(inm,"%d %d %d %d %d %d %d %d %d\n",hexTab[9*k],hexTab[9*k+1],hexTab[9*k+2],
//            hexTab[9*k+3],hexTab[9*k+4],hexTab[9*k+5],hexTab[9*k+6],hexTab[9*k+7],hexTab[9*k+8] );
//  }
//
//  fprintf(inm,"\n\nQuadrilaterals %d\n",quadNumber);
//  for (int k=1; k<=mmgMesh->nquad; k++) {
//    MMG5_pQuad pq = &mmgMesh->quadra[k];
//    fprintf(inm,"%d %d %d %d %d\n",pq->v[0],pq->v[1],pq->v[2],pq->v[3],pq->ref);
//  }
//
//
//  fprintf(inm,"\n\nEnd\n");
//
//  fclose(inm);

  /** Step 5: converts hexa into a MMG5 tetrahedral mesh */
  ier = H2T_libhex2tet(mmgMesh,hexTab,nbHex);

  if ( ier != H2T_STRONGFAILURE ) {
    MMG3D_saveMesh(mmgMesh,fileout);
  }
  else {
    printf("Hex2tet Fail: unable to save mesh.\n");
  }

  free(hexTab);

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return ier;
}
