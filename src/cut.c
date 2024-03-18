/**
 * \file cut.c
 * \brief Functions to cut the hexa mesh into a tetra one
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
#include <memory.h>

#define H2T_MAXTET_ERROR_MESSAGE(func,line,nemax,ncut,nhex) do          \
  {                                                                     \
    fprintf(stdout,"%s:%d: max number of tet reached (%" MMG5_PRId "). %d/%d hexa treated.\n", \
            (func),(line),(nemax),(ncut),(nhex));                       \
  } while(0)


/**
 * \param mesh pointer toward the tetrahedral mesh
 * \param hed edge hash table
 * \param v0 hexa vertex
 * \param v1 hexa vertex
 * \param v2 hexa vertex
 * \param v3 hexa vertex
 * \param ncut pointer toward the number of hexa cutted (to update)
 * \param nhex total number of hexa
 *
 * \return 0 if fail, 1 if success
 *
 * Add a tetrahedron of vertices p[0]\@4 and print error messages if fail.
 *
 */
static
int H2T_Add_tetra(MMG5_pMesh mesh,int v0,int v1,int v2,int v3,int ref,int ncut,int nhex) {
  int iel;

  iel =  MMG3D_Add_tetrahedron ( mesh,v0,v1,v2,v3,ref);

  if ( !iel ) {
    H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut,nhex);
    return 0;
  }
  if ( iel < 0 ) {
    printf("  ## Error:%s: tetra %d %d %d %d: wrong orientation\n",
           __func__,v0,v1,v2,v3);
    return 0;
  }
  //mesh->tetra[iel].mark=0;

  return 1;
}

/**
 * \param mesh pointer toward the tetrahedral mesh
 * \param hed edge hash table
 * \param p hexa vertices
 * \param ncut pointer toward the number of hexa cutted (to update)
 * \param nhex total number of hexa
 *
 * \return 0 if fail, 1 if success
 *
 * Cut the hexahedron of vertices p[0]\@7 into 6 tetrahedra.
 *
 */
static
int H2T_decouphex(MMG5_pMesh mesh, pHedge hed,int* p,int ref,int *ncut,int nhex) {

  /** Creation of the first tetra (0,1,3,7) */
  if ( !H2T_Add_tetra(mesh,p[0],p[1],p[3],p[7],ref,*ncut,nhex) ) return 0;

  /** Creation of the second tetra (7,2,6,1) */
  if ( !H2T_Add_tetra(mesh,p[7],p[2],p[6],p[1],ref,*ncut,nhex) ) return 0;

  /** Creation of the third tetra (1,4,5,7) */
  if ( !H2T_Add_tetra(mesh,p[1],p[4],p[5],p[7],ref,*ncut,nhex) ) return 0;

  /** Creation of the fourth tetra (7,4,0,1) */
  if ( !H2T_Add_tetra(mesh,p[7],p[4],p[0],p[1],ref,*ncut,nhex) ) return 0;

  /** Creation of the fourth tetra (1,6,7,5) */
  if ( !H2T_Add_tetra(mesh,p[1],p[6],p[7],p[5],ref,*ncut,nhex) ) return 0;

  /** Creation of the sixth tetra (1,3,7,2) */
  if ( !H2T_Add_tetra(mesh,p[1],p[3],p[7],p[2],ref,*ncut,nhex) ) return 0;

  /** Add edges to the hashtable */
  if ( !H2T_edgePut(hed,p[0],p[7],2) ) return 0;
  if ( !H2T_edgePut(hed,p[1],p[3],2) ) return 0;
  if ( !H2T_edgePut(hed,p[2],p[7],2) ) return 0;
  if ( !H2T_edgePut(hed,p[1],p[6],2) ) return 0;
  if ( !H2T_edgePut(hed,p[1],p[4],2) ) return 0;
  if ( !H2T_edgePut(hed,p[5],p[7],2) ) return 0;

  ++(*ncut);

  return 1;
}

/**
 * \param ph hexa vertices
 * \param nu1 First extremity of 1 created diagonale inside the hexa
 * \param nu2 Second extremity of 1 created diagonale inside the hexa
 * \param hed edge hashtable
 *
 * \return 1 if diagonal with nu1 extremity other than nu1-nu2 exists
 *
 * Check if we have already hashed at least 1 diagonal issue from nu1
 * (of the hexa faces) without taking into account the nu1-nu2
 * diagonal (default case).
 *
 */
static int H2T_checkcase(int ph[8],int nu1,int nu2,pHedge hed) {
  int i,nu3;

  for ( i=0; i<3; i++ ) {
    /** Check the existence of edge nu1-nu3 with nu3 the point
     * opposite to nu1 inside the hexa faces */
    nu3 = H2T_hied[nu1][i];

    /** Do not check the edge from which we come */
    if ( nu3==nu2 ) continue;

    /** Test the edge existence */
    if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
  }

  if ( i<3 ) return 1;
  else return 0;
}

/**
 * \param ph hexa vertices
 * \param nu1 First extremity of a non created inside the hexa
 * \param nu2 Second extremity of a non created diagonale inside the hexa
 * \param hed edge hashtable
 *
 * \return 1 if at least 2 edges of the opposite case to the default case exist.
 *
 * Check if we have already hashed at least 1 diagonal issue from nu1
 * (of the hexa faces) without taking into account the nu1-nu2
 * diagonal.
 *
 */
static int H2T_checkcaseopp(int ph[8],int nu1,int nu2,pHedge hed) {
  int i,nu3,nu4;

  /** Point opposite to nu1 in the hexa */
  nu4 = H2T_hop[nu1];

  for ( i=0; i<3; i++ ) {
    /** Check the existence of edge nu4-nu3 with nu3 the point
     * opposite to nu1 inside the hexa faces */
    nu3 = H2T_hied[nu4][i];

    /** Do not check the edge hop[nu1]-hop[nu2] */
    if ( nu3==H2T_hop[nu2] ) continue;

    /** Test the edge existence */
    if ( H2T_edgePoint(hed,ph[nu4],ph[nu3]) ) break;
  }

  if ( i<3 ) return 1;
  else return 0;
}

/**
 *
 * Check if the adja array allow to reach all the hexa by adjacency
 * form the first one.
 *
 */
int H2T_chkAdja(MMG5_pMesh mesh,int* listhexa,MMG5_int* adjahex,int nhex) {
  int *list,*mark;
  int icurc,ipil,iadr,adj,count,k,i;

  /** Allocs */
  list = mark = NULL;
  list = (int*) calloc(7*nhex+1,sizeof(long long));
  assert(list);
  mark = (int*) calloc(nhex+1,sizeof(long long));
  assert(mark);


  /** Stack initialization: mark the frist hexa as seen and append its
   * adjacents to the stack */
  icurc   = 0;
  ipil    = 0;
  mark[1] = 1;
  iadr    = 1;

  for ( i=0; i<6; i++ ) {
    adj  = adjahex[iadr + i];

    if ( !adj ) continue;

    list[ipil++] = adj;
  }

  /** Process the next hexa of the stack: mark it as seen and append its adjacents */
  while ( icurc++ < ipil ) {
    k = list[icurc-1]/6;

    if ( !k ) continue;

    mark[k] = 1;

    for ( i=0; i<6; i++ ) {
      iadr = 6*(k-1)+1;
      adj = adjahex[iadr + i];

      if ( !adj ) continue;

      /* Add only once each tetra */
      if ( mark[adj/6] )
        continue;
      else {
        mark[adj/6] = -1;
      }

      list[ipil++] = adj;
    }
  }

  /** Count the number of marked tetra (already seen) */
  count = 0;
  for ( i=1; i<=nhex; ++i ) {
    if ( mark[i] > 0 ) ++count;
  }

  /** Free the memory */
  free(list); list = NULL;
  free(mark); mark = NULL;

  return count;
}

/*
 * Hexahedron
 *
 *
 *   4----------5          .----------.           .----------.
 *   |\         |\         |\         |\          |\      _/ |\
 *   | \        | \        | \     3  | \         | \   _/   | \
 *   |  \       |  \       |  \  5    |  \        |  \ /     |  \
 *   |   7------+---6      |   .------+---.       |  /.------+---.
 *   |   |      |   |      | 1 |      | 4 |       | / |\_    |   |
 *   0---+------1   |      .---+------.   |       ./--+--\_--.   |
 *    \  |       \  |       \  |     2 \  |        \ \|__  \_ \  |
 *     \ |        \ |        \ |  0     \ |         \ |  \__ \_\ |
 *      \|         \|         \|         \|          \|     \__ \|
 *       3----------2          .----------.           .----------.
 *
 */
int H2T_cuthex(MMG5_pMesh mesh,pHedge hed,int* listhexa,MMG5_int* adjahex,int nhex) {
  MMG5_pPoint    ppt;
  int            i,ih,nu1,nu2,nu3,nu4,adj,icas0,icasopp,nncut;
  int            *list,*mark,p[8],ipil,icurc,iface,iadr;
  int            iel,ip,ph[8];
  double         c[3];
  int            ddebug,ncut;
  MMG5_int k;

  if ( mesh->info.ddebug ) {
    int count = H2T_chkAdja(mesh,listhexa,adjahex,nhex);
    printf("Number of hexa reached by adjacency: %d/%d\n",count,nhex);
  }

  /* Alloc */
  list = (int*) calloc(7*nhex+1,sizeof(long long));
  assert(list);
  mark = (int*) calloc(nhex+1,sizeof(long long));
  assert(mark);

  /* Stack initialization */
  mark[1] = -1;

  for ( ih=0; ih<8; ih++ ) p[ih] = listhexa[9+ih];

  if ( mesh->info.ddebug )
    printf(" First hexa %d %d %d %d %d %d %d %d\n",p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);

  ncut = 0;
  if ( !H2T_decouphex(mesh,hed,p,listhexa[9+8],&ncut,nhex) ) return 0;

  icurc = 0;
  ipil  = 0;
  for ( i=0; i<6; i++ ) {
    iadr = 1;
    adj = adjahex[iadr + i];
    if ( !adj) continue;
    list[ipil++] = adj;
  }

  while(icurc++ < ipil) {
    k = list[icurc-1]/6;

    if ( !k ) continue;

    /** Store the hexa */
    for(ih=0;ih<8;ih++) ph[ih] = listhexa[9*k+ih];
    if ( mesh->info.ddebug ) {
      printf("ph %d %d %d %d %d %d %d %d\n",
             ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
    }

    mark[k] = -1;

    iface = list[icurc-1]%6;
    if ( iface < 0 ) {
      printf(" ## Error: %s: wrong adjacency (%d)\n",__func__,iface);
      return 0;
    }
    ddebug=0;

    nu1 = H2T_hidir[iface][0];
    nu2 = H2T_hidir[iface][2];

    /** Test if the diagonal 0-2 of iface exist */
    if ( H2T_edgePoint(hed,ph[nu1],ph[nu2]) ) {

      /** Test if the opposite diagonal on the opposite face exist */
      nu3 = H2T_hidirop[iface][0];
      nu4 = H2T_hidirop[iface][2];
      if ( H2T_edgePoint(hed,ph[nu3],ph[nu4]) ) {
        /** If yes, mark the hexa */
        mark[k] = -10;
        continue;
      }
      if ( iface==1 || iface==5 ) {
        //find if other edge with ph->v[MMG_hidir[iface][0]], if yes->renum
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if (  icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //debug check
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) break;
          }
          assert(i==3);
          //printf("iface %d: another edge founded---> renum\n",iface);
          if ( iface==1 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          }
          if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0;
        } else {
	  if ( !H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex) ) return 0;
        }
      } else if ( iface==4 ) {
        icas0 = H2T_checkcase(ph,nu2,nu1,hed);
        icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu1,nu2,hed);
          icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu1][i];
            if ( nu3==nu2 ) continue;
            if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
          }
          assert(i==3);
          //printf("iface 4 another edge founded ---> renum\n");
          p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
          p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0;
        } else {
          if ( !H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex) ) return 0;
        }
      } else {
        if ( ddebug)  printf("face %d renumbering\n",iface);//iface 0,2,3
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) break;
          }
          assert(i==3);
          icas0=1;
        }
        if ( ddebug) printf("icas %d\n",icas0);
        switch(iface) {
        case(0):
          if ( icas0 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          }
          break;
        case(2):
          if ( icas0 ) {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          break;
        case(3):
          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          }
          break;
        }

        if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0; 
      }
    }  else if ( H2T_edgePoint(hed,ph[H2T_hidir[iface][1]],ph[H2T_hidir[iface][3]]) ) {
      /** The edge 1-3 (opposite to 0-2) of iface exist */
      nu1 = H2T_hidir[iface][1];
      nu2 = H2T_hidir[iface][3];

      nu3 = H2T_hidirop[iface][1];
      nu4 = H2T_hidirop[iface][3];
      /** Test if the opposite edge on the opposite face exist */
      if ( H2T_edgePoint(hed,ph[nu3],ph[nu4]) ) {
        /** nothing to do */
        mark[k] = -10;
        continue;
      }
      if ( iface==0 || iface==3 ) {
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( ddebug) printf("icas0 %d\n",icas0);
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu2][i];
            if ( nu3==nu1 ) continue;
            if ( H2T_edgePoint(hed,ph[nu2],ph[nu3]) ) {
              break;
            }
          }
          assert(i==3);
          if ( iface==0 ) {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0;
        } else {
          if ( !H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex) ) return 0;
        }
      } else if ( iface==2 ) {
        icas0 = H2T_checkcase(ph,nu2,nu1,hed);
        icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu1,nu2,hed);
          icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        if ( icas0 ) {
          //check debug
          for ( i=0; i<3; i++ ) {
            nu3 = H2T_hied[nu1][i];
            if ( nu3==nu2 ) continue;
            if ( H2T_edgePoint(hed,ph[nu1],ph[nu3]) ) break;
          }
          assert(i==3);
          p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
          p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0;
        } else {
          if ( !H2T_decouphex(mesh,hed,ph,listhexa[9*k+8],&ncut,nhex) ) return 0;
        }
      }
      else {
        if ( ddebug)  printf("face %d renumbering\n",iface);//iface 1,4,5
        icas0 = H2T_checkcase(ph,nu1,nu2,hed);
        icasopp = H2T_checkcaseopp(ph,nu1,nu2,hed);
        if ( icas0 || icasopp ) {
          icas0 = H2T_checkcase(ph,nu2,nu1,hed);
          icasopp = H2T_checkcaseopp(ph,nu2,nu1,hed);
          if ( icas0 || icasopp ) {
            mark[k] = -10;
            continue;
          }
          icas0 = 1;
        }
        switch(iface) {
        case(1):
          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          break;
        case(4):
          if ( ddebug ) {
            printf("at the beginning %" MMG5_PRId" : %d %d %d %d %d %d %d %d\n",
                   k,ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
          }

          if ( icas0 ) {
            p[0] = ph[1]; p[1] = ph[2]; p[2] = ph[3]; p[3] = ph[0];
            p[4] = ph[5]; p[5] = ph[6]; p[6] = ph[7]; p[7] = ph[4];
          } else {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          }
          if ( ddebug ) {
            printf("at the end %" MMG5_PRId" : %d %d %d %d %d %d %d %d\n",
                   k,p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7]);
          }
          break;
        case(5):
          if ( icas0 ) {
            p[0] = ph[2]; p[1] = ph[3]; p[2] = ph[0]; p[3] = ph[1];
            p[4] = ph[6]; p[5] = ph[7]; p[6] = ph[4]; p[7] = ph[5];
          } else {
            p[0] = ph[3]; p[1] = ph[0]; p[2] = ph[1]; p[3] = ph[2];
            p[4] = ph[7]; p[5] = ph[4]; p[6] = ph[5]; p[7] = ph[6];
          }
          break;
        }
        if ( !H2T_decouphex(mesh,hed,p,listhexa[9*k+8],&ncut,nhex) ) return 0;
      }
    } else {
      /** We have no edge on iface: not normal because the adjacent
       * through iface has been cutted */
      printf("Error: %s: no created diagonal on face %d of hexa: %d %d and %d %d tested\n",
             __func__,iface,ph[H2T_hidir[iface][1]],ph[H2T_hidir[iface][3]],
             ph[H2T_hidir[iface][0]],ph[H2T_hidir[iface][2]]);
      printf("hexa: %d %d %d %d %d %d %d %d\n",ph[0],ph[1],ph[2],ph[3],ph[4],ph[5],ph[6],ph[7]);
      return 0;
    }

    /** Append the adjacents of the current tetra to the stack */
    for ( i=0; i<6; i++ ) {
      iadr = 6*(k-1)+1;
      adj = adjahex[iadr + i];

      if ( !adj ) continue;

      if ( mark[adj/6] )
        continue;
      else
        mark[adj/6] = 10;

      list[ipil++] = adj;

      if ( mesh->info.ddebug )
        printf("k=%" MMG5_PRId": stack append: hexa %d (iface %d) in %d -- through face %d\n",
               k,adj/6,adj%6,ipil-1,i);
    }
  }

  if ( ncut < nhex ) {
    printf("     %d/%d Hexa succesfully cutted\n",ncut,nhex);
    printf("--> Try to cut the remaining hexa using point insertion\n");
  }

  /** Treat the remaining hexa */
  nncut = 0;
  for ( k=1; k<=nhex; k++ ) {
    if ( mark[k]==-1 ) continue;
    for(i=0;i<8;i++) ph[i] = listhexa[9*k+i];
    nncut++;
    /** create new vertex */
    c[0] = c[1] = c[2] = 0.;
    for ( i=0; i<8; i++ ) {
      ppt = &mesh->point[ph[i]];
      c[0] += ppt->c[0];        c[1] += ppt->c[1]; c[2] += ppt->c[2];
    }
    c[0] /= 8.; c[1] /= 8.; c[2] /= 8.;

    ip = MMG3D_Add_vertex(mesh,c[0],c[1],c[2],0);
    if ( !ip ) return 0;

    /** create 2 tets per faces */
    int ref = listhexa[9*k+8];

    for ( i=0 ;i<6; i++ ) {
      nu1 = H2T_hidir[i][0];
      nu2 = H2T_hidir[i][2];
      if ( H2T_edgePoint(hed,ph[nu1],ph[nu2]) ) {
        if ( mesh->ne+1 >= mesh->nemax ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][1]] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][1]],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[H2T_hidir[i][3]],ph[nu2] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[H2T_hidir[i][3]],ph[nu2],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }
      } else {
        nu1 = H2T_hidir[i][1];
        nu2 = H2T_hidir[i][3];
        if ( !H2T_edgePoint(hed,ph[nu1],ph[nu2])) { 
	  iel = H2T_edgePut(hed,ph[nu1],ph[nu2],2);
	}
	if ( !iel ) return 0;
        if ( mesh->ne+1 >= mesh->nemax ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[H2T_hidir[i][0]],ph[nu2] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[H2T_hidir[i][0]],ph[nu2],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }

        /* Tetra ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][2]] */
        iel =  MMG3D_Add_tetrahedron ( mesh,ip,ph[nu1],ph[nu2],ph[H2T_hidir[i][2]],ref );
        if ( !iel ) {
          H2T_MAXTET_ERROR_MESSAGE(__func__,__LINE__,mesh->nemax,ncut+nncut,nhex);
          fprintf(stdout,"%d new points.\n",nncut);
          return 0;
        }
      }
    }
  }
  if ( nncut) fprintf(stdout,"  $$ %8d ADDED VERTEX\n",nncut);
  return 1;
}
