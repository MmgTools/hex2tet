/**
 * \file hashhexa.c
 * \brief Hashing functions.
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
#include <limits.h>

/** Constant for hash table key computation */
#define H2T_KA 31
/** Constant for hash table key computation */
#define H2T_KB 57
/** Constant for hash table key computation */
#define H2T_KC 79

/** Find the minimum value and location of v0, v1, v2 and v3 */
#define H2T_MINVAL_AND_LOC(v0,v1,v2,v3,minval,minloc) do \
  {                                                      \
                                                         \
    minval = v0; minloc = 0;                             \
    if ( v1 < minval ) {                                 \
      minval = v1; minloc = 1;                           \
    }                                                    \
    if ( v2 < minval ) {                                 \
      minval = v2; minloc = 2;                           \
    }                                                    \
    if ( v3 < minval ) {                                 \
      minval = v3; minloc = 3;                           \
    }                                                    \
  } while (0)

/**
 * \param listhexa pointer toward the list of hexa (hexa number i is stored in listhexa[9*i])
 * \param adjhex adjacency array to fill
 * \param nhex number of hexa in listhexa
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Create table of adjacency.
 *
 */
int H2T_hashHexa(int* listhexa,MMG5_int* adjahex,int nhex) {
  MMG5_int         k,kk,pp,l,ll,mins,mins1,opps,opps1,sum,sum1,iadr,*link;
  int              *hcode,hsize,imins,imins1;
  int              ph[8],ph1[8];
  long int         inival;
  unsigned char    i,ii,iii,i1,i2,i3,i4;
  unsigned int     key;

  /* default */
  fprintf(stdout,"  ** SETTING HEXA ADJACENCIES\n\n");
  fflush(stdout);

  /* memory alloc */
  hcode = (int*)calloc(nhex+1,sizeof(long long));
  assert ( hcode );
  link  = adjahex;
  hsize = nhex;

  /* init */
  inival = INT_MAX;
  for (k=0; k<=nhex; k++)
    hcode[k] = -inival;

  /* build hash table */
  for ( k=1; k<=nhex; k++ ) {
    iadr = 9*k;
    for ( iii=0; iii<8; iii++ )
      ph[iii] = listhexa[iadr+iii];

    for ( i=0; i<6; i++ ) {
      i1 = H2T_hidir[i][0];
      i2 = H2T_hidir[i][1];
      i3 = H2T_hidir[i][2];
      i4 = H2T_hidir[i][3];

      H2T_MINVAL_AND_LOC(ph[i1],ph[i2],ph[i3],ph[i4],mins,imins);
      opps = ph[H2T_hopp[i][imins]];

      /* compute key */
      sum = ph[i1] + ph[i2] + ph[i3] + ph[i4];
      key = H2T_KA*mins + H2T_KB*opps + H2T_KC*sum;
      key = key % hsize + 1;

      /* insert */
      iadr = 6*(k-1) + i+1;
      link[iadr] = hcode[key];
      hcode[key] = -iadr;
    }
  }

  /* set adjacency */
  for (l=6*nhex; l>0; l--) {

    if ( link[l] >= 0 )  continue;

    /* current element */
    k = ((l-1) /6) + 1;
    i = (l-1) % 6;
    i1 = H2T_hidir[i][0];
    i2 = H2T_hidir[i][1];
    i3 = H2T_hidir[i][2];
    i4 = H2T_hidir[i][3];
    iadr = 9*k;
    for(iii=0 ; iii<8 ; iii++)
      ph[iii] = listhexa[iadr+iii];
    sum  = ph[i1] + ph[i2] + ph[i3] + ph[i4];

    H2T_MINVAL_AND_LOC(ph[i1],ph[i2],ph[i3],ph[i4],mins,imins);
    opps = ph[H2T_hopp[i][imins]];

    /* accross link */
    ll = -link[l];
    pp = 0;
    link[l] = 0;
    while ( ll != inival ) {
      kk = ((ll-1) /6) + 1;
      ii = (ll-1) % 6;
      i1 = H2T_hidir[ii][0];
      i2 = H2T_hidir[ii][1];
      i3 = H2T_hidir[ii][2];
      i4 = H2T_hidir[ii][3];
      iadr = 9*kk;
      for(iii=0 ; iii<8 ; iii++)
        ph1[iii] = listhexa[iadr+iii];
      sum1 = ph1[i1] + ph1[i2] + ph1[i3] + ph1[i4];
      if ( sum1 == sum ) {

        H2T_MINVAL_AND_LOC(ph1[i1],ph1[i2],ph1[i3],ph1[i4],mins1,imins1);
        opps1 = ph1[H2T_hopp[ii][imins1]];

        if ( mins1 == mins ) {
          if ( opps1 == opps ) {
            /* adjacent found */
            if ( pp != 0 )  link[pp] = link[ll];
            link[l] = 6*kk + ii;
            link[ll]= 6*k + i;
            break;
          }
        }
      }
      pp = ll;
      ll = -link[ll];
    }
  }

  free(hcode); hcode=NULL;

  return 1;
}

/**
 * \param hash edge hash table
 * \param a edge extremity
 * \param b edge extremity
 *
 * \return the edge flag if the edge exists, 0 otherwise.
 *
 * Find the flag of the edge a-b of the hash table.
 *
 */
int H2T_edgePoint(pHedge hash,int a,int b) {
  MMG5_int  key,mins,maxs;
  hedge     *ha;

  /* compute key */
  mins = a;
  maxs = b;
  if ( a > b ) {
    mins = b;
    maxs = a;
  }
  key = H2T_KA*mins + H2T_KB*maxs;
  key = key % hash->size;
  ha  = &hash->item[key];

  if ( !ha->min )  return 0;
  else if ( ha->min == mins && ha->max == maxs ) {
    return ha->iel;
  }
  else if ( ha->nxt ) {
    do {
      ha = &hash->item[ha->nxt];
      if ( ha->min == mins && ha->max == maxs ) {
        return ha->iel;
      }
    }
    while ( ha->nxt && ha->nxt < hash->nhmax );
  }
  return 0;
}

/**
 * \param hash edge hash table
 * \param a edge extremity
 * \param b edge extremity
 * \param np flag to set on edge
 *
 * \return the edge flag if the edge exists in the hash table,1 if we
 * add the edge, 0 if fail
 *
 * Hash the edge a-b and flag it with the integer \a np.
 *
 */
int H2T_edgePut(pHedge hash,int a,int b,int np) {
  MMG5_int  key,mins,maxs;
  hedge     *ha;

  mins = a;
  maxs = b;
  if ( a > b ) {
    mins = b;
    maxs = a;
  }
  key = H2T_KA*mins + H2T_KB*maxs;
  key = key % hash->size;
  ha  = &hash->item[key];

  if ( ha->min ) {
    /* Same edge */
    if ( ha->min == mins && ha->max == maxs ) {
      return ha->iel;
    }
    else {
      while ( ha->nxt && ha->nxt < hash->nhmax ) {
        ha = &hash->item[ha->nxt];
        if ( ha->min == mins && ha->max == maxs ) {
          return ha->iel;
        }
      }
    }
    ha->nxt = hash->hnxt;
    ha      = &hash->item[hash->hnxt];
    ++hash->hnxt;
    if ( hash->hnxt >= hash->nhmax ) {
      fprintf(stdout,"  ## Memory alloc problem (edge): %" MMG5_PRId "\n",hash->nhmax);
      return 0;
    }
  }

  /* insert */
  ha->min = mins;
  ha->max = maxs;
  ha->iel = np;
  ha->nxt = 0;

  return 1;
}
