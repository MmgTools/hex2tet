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

/** Constant for hash table key computation (as in Mmg) */
#define H2T_MMG5_KA 7
/** Constant for hash table key computation (as in Mmg) */
#define H2T_MMG5_KB 11
/** Constant for table reallocation */
#define H2T_GAP 1.2

/* Macros */
#define H2T_MAX(a,b) (((a) > (b)) ? (a) : (b))
#define H2T_MIN(a,b) (((a) < (b)) ? (a) : (b))


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
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table of edges.
 * \param hsiz initial size of hash table.
 * \param hmax maximal size of hash table.
 * \return 1 if success, 0 if fail.
 *
 * Hash edges or faces.
 *
 * \remark Modeled after MMG5_hashNew.
 */
int H2T_hashNew(MMG5_pMesh mesh,MMG5_Hash *hash,int hsiz,int hmax) {
  int   k;

  /* adjust hash table params */
  hash->siz  = hsiz+1;
  hash->max  = hmax + 2;
  hash->nxt  = hash->siz;

  hash->item = (MMG5_hedge*)calloc(hash->max+1,sizeof(MMG5_hedge));

  for (k=hash->siz; k<hash->max; k++)
    hash->item[k].nxt = k+1;

  return 1;
}

/**
 * \param mesh pointer toward the mesh.
 * \param hash pointer toward the hash table to fill.
 * \param ia first vertex of face to hash.
 * \param ib second vertex of face to hash.
 * \param ic third vertex of face to hash.
 * \param k index of face to hash.
 *
 * \return 0 if fail, -1 if the face is newly hashed, index of the first face
 * hashed if another face with same vertices exist.
 *
 * \remark Modeled after MMG5_hashFace.
 *
 **/
int H2T_hashFace(MMG5_pMesh mesh,MMG5_Hash *hash,int ia,int ib,int ic,int k) {
  MMG5_hedge     *ph;
  int        key,mins,maxs,sum,j;

  mins = H2T_MIN(ia,H2T_MIN(ib,ic));
  maxs = H2T_MAX(ia,H2T_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (H2T_MMG5_KA*mins + H2T_MMG5_KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return ph->k;
    else {
      while ( ph->nxt && ph->nxt < hash->max ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return ph->k;
      }
    }
    ph->nxt = hash->nxt;
    ph      = &hash->item[hash->nxt];
    ph->a   = mins;
    ph->b   = maxs;
    ph->s   = sum;
    ph->k   = k;
    hash->nxt = ph->nxt;
    ph->nxt = 0;

#warning no hash table reallocation
    if ( hash->nxt >= hash->max ) {
      hash->max *= H2T_GAP;
      hash->item = (MMG5_hedge*)realloc(hash->item,hash->max*sizeof(MMG5_hedge));
      for (j=hash->nxt; j<hash->max; j++) {
        ph = &hash->item[j];
        ph->nxt = j+1;
        ph->a = ph->b = ph->s = ph->k = 0;
      }
    }
    return -1;
  }

  /* insert new face */
  ph->a = mins;
  ph->b = maxs;
  ph->s = sum;
  ph->k = k;
  ph->nxt = 0;

  return -1;
}

/**.
 * \param hash pointer toward the hash table.
 * \param ia first triangle vertex
 * \param ib second triangle vertex
 * \param ic third triangle vertex
 *
 * \return 0 if failed, the index of the triangle otherwise.
 *
 * Look for the triangle with vertices (ia,ib,ic) in the hash table, and return
 * its index.
 *
 * \remark Modeled after MMG5_hashGetFace.
 *
 */
int H2T_hashGetFace(MMG5_Hash *hash,int ia,int ib,int ic) {
  MMG5_hedge  *ph;
  int     key,mins,maxs,sum;

  if ( !hash->item )  return 0;

  mins = H2T_MIN(ia,H2T_MIN(ib,ic));
  maxs = H2T_MAX(ia,H2T_MAX(ib,ic));

  /* compute key */
  sum = ia + ib + ic;
  key = (H2T_MMG5_KA*mins + H2T_MMG5_KB*maxs) % hash->siz;
  ph  = &hash->item[key];

  if ( ph->a ) {
    if ( ph->a == mins && ph->b == maxs && ph->s == sum )
      return ph->k;
    else {
      while ( ph->nxt ) {
        ph = &hash->item[ph->nxt];
        if ( ph->a == mins && ph->b == maxs && ph->s == sum )  return ph->k;
      }
    }
  }

  return 0;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * Hash all possible triangles that could be created from input boundary
 * quadrilaterals (4*nquad) and store them with index 4*quad_id+i.
 *
 */
int H2T_hashQuad(MMG5_pMesh mesh,MMG5_Hash *hash) {
  MMG5_pQuad pq;
  int k;

  hash->item = NULL;

  if ( ! H2T_hashNew(mesh,hash,4*mesh->nquad,1.51*4*mesh->nquad) ) return 0;

  int nh = 0;
  for( k = 1; k <= mesh->nquad; k++) {
    pq = &mesh->quadra[k];

    /* Hash all the four possible triangle cut */
    if ( !H2T_hashFace(mesh,hash,pq->v[0],pq->v[1],pq->v[2],4*k) ) {
      free(hash->item);
      return 0;
    }
    ++nh;
    if ( !H2T_hashFace(mesh,hash,pq->v[2],pq->v[3],pq->v[0],4*k+1) ) {
      free(hash->item);
      return 0;
    }
    ++nh;
    if ( !H2T_hashFace(mesh,hash,pq->v[1],pq->v[2],pq->v[3],4*k+2) ) {
      free(hash->item);
      return 0;
    }
    ++nh;
    if ( !H2T_hashFace(mesh,hash,pq->v[3],pq->v[0],pq->v[1],4*k+3) ) {
      free(hash->item);
      return 0;
    }
    ++nh;
  }
  /* Check that all the possible triangles have been hashed */
  assert(nh == 2*mesh->nt);

  return 1;
}

/**
 * \param mesh pointer toward the mesh structure.
 * \param hash pointer toward the hash table.
 *
 * \return 0 if failed, 1 otherwise.
 *
 * If a boundary triangles is found in the hash table, add it to the mesh and
 * give it the reference of the corresponding input boundary quadrilateral.
 *
 */
int H2T_hashGetRef(MMG5_pMesh mesh,MMG5_Hash *hash) {
  MMG5_pTetra pt;
  MMG5_int         *adja;
  int         ie,i,nt,ia,ib,ic,k,ref;

  /* Create tetrahedra adjacency */
  if( !MMG3D_hashTetra(mesh,0) ) return 0;
  assert(mesh->adja);

  /* Look for boundary triangles in the hash table */
  nt = 0;
  for( ie = 1; ie <= mesh->ne; ie++ ){
    pt = &mesh->tetra[ie];
    adja = &mesh->adja[4*(ie-1)+1];
    for( i = 0; i < 4; i++ ) {
      if( adja[i] ) continue;
      ia = pt->v[H2T_tidir[i][0]];
      ib = pt->v[H2T_tidir[i][1]];
      ic = pt->v[H2T_tidir[i][2]];
      assert( ia && ib && ic );
      k = H2T_hashGetFace(hash,ia,ib,ic);
      if( !k ) continue; /* face not found */
      k /= 4;
      ref = mesh->quadra[k].ref;
      if( !MMG3D_Set_triangle(mesh,ia,ib,ic,ref,++nt) )
        return 0;
    }
  }
  /* Check that you found all the triangles coming from a cut quadrilateral */
  assert(nt == mesh->nt);

  /* Free hash table */
  free(hash->item);

#warning quads will be freed by MMG3D_Free_all
  /* Clean quadrilaterals */
  mesh->nquad = 0;

  return 1;
}

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
