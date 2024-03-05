/**
 * \file chk.c
 * \brief Hexa volume and orientation computations.
 *
 * \author Cecile Dobrzynski (Bx INP/Inria/UBordeaux)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

/**
 * \param mesh pointer toward the mesh structure.
 * \param iel element index
 *
 * \return the voume of the tetra.
 *
 * Compute the tetra volume from the tetra index.
 *
 */
double H2T_voltet(MMG5_pMesh mesh,int iel) {
  MMG5_pTetra pt;

  pt = &mesh->tetra[iel];

  return H2T_quickvol(mesh->point[pt->v[0]].c,mesh->point[pt->v[1]].c,
                      mesh->point[pt->v[2]].c,mesh->point[pt->v[3]].c);
}

/**
 * \param c1 coor of the first vertex
 * \param c2 coor of the second vertex
 * \param c3 coor of the third vertex
 * \param c4 coor of the fourth vertex
 *
 * \return the voume of the tetra
 *
 * Compute the tetra volume from its vertices coor.
 *
 */
double H2T_quickvol(double *c1,double *c2,double *c3,double *c4) {
  double   ax,ay,az,bx,by,bz,vol;

  ax = c3[0] - c1[0];
  ay = c3[1] - c1[1];
  az = c3[2] - c1[2];

  bx = c4[0] - c1[0];
  by = c4[1] - c1[1];
  bz = c4[2] - c1[2];

  vol = (c2[0]-c1[0]) * (ay*bz - az*by) \
    + (c2[1]-c1[1]) * (az*bx - ax*bz) \
    + (c2[2]-c1[2]) * (ax*by - ay*bx);

  return vol;
}

/**
 * \param mmgMesh pointer toward the mesh structure.
 * \param hexa hexa array
 * \param nhex number of hexa in the hexa array.
 *
 * \return the number of reorientated hexa
 *
 * Count the number of hexa with bad orientation and permute their
 * vertices in order to have only "good" orientated hexa.
 *
 */
int H2T_chkorient(MMG5_pMesh mmgMesh,int* hexa,int nhex) {
  int     nbado,k,i,ph[8];
  double  volref,volhex;
  MMG5_int iadr;

  nbado = 0;
  volref = 1;

  for (k=1; k<=nhex; k++) {
    iadr = (MMG5_int)9*k;
    for(i=0 ; i<8 ; i++)
      ph[i] = hexa[iadr+i];

    /** check the orientability of the hexahedra : vol of tet p0 p1 p3 p4 */
    volhex = H2T_quickvol(mmgMesh->point[ph[0]].c,mmgMesh->point[ph[1]].c
                          ,mmgMesh->point[ph[3]].c,mmgMesh->point[ph[4]].c);

    if ( volref*volhex < 0 ) {
      nbado++;
      hexa[iadr + 3] = ph[1];
      hexa[iadr + 1] = ph[3];
      hexa[iadr + 5] = ph[7];
      hexa[iadr + 7] = ph[5];
    }
  }
  return nbado;
}
