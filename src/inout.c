/**
 * \brief I/O functions
 *
 * \author Cecile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (Inria)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define H2T_HEX_LINE_FEED 0x0a
#define H2T_HEX_APOS 0x27
#define H2T_HEX_LEFTPAR 0x28
#define H2T_HEX_RIGHTPAR 0x29
#define H2T_HEX_COMMA 0x2c
#define H2T_HEX_ZERO 0x30
#define H2T_HEX_NINE 0x39

static inline
MMG5_int H2T_npy_point_index(int i,int j,int k, int *n) {
  return (n[1]+1)*(n[2]+1)*i + (n[2]+1)*j + k + 1;
}

int H2T_loadNpy(MMG5_pMesh mmgMesh, int** tabhex, char* filename) {

  FILE* inm;
  char buffer = 0;
  char* str = NULL;
  size_t dataSize = 0;
  int pos1, pos2, dim = 0, t[3];
  MMG5_int nhex, np, ne, i, j, k, ref, pos;

  /* Input data and creation of the hexa array */
  if( !(inm = fopen(mmgMesh->namein,"rb")) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",mmgMesh->namein);
    return 0;
  }

  fprintf(stdout,"  %%%% %s OPENED\n",mmgMesh->namein);

  /* Reach beginning of dimension specification in header */
  while (!(buffer == H2T_HEX_LEFTPAR)) {
    fread(&buffer,sizeof(buffer),1,inm);
  }
  
  /* Read array sizes */
  while (!(buffer == H2T_HEX_RIGHTPAR)) {
  
    pos1 = ftell(inm);
    do {
      fread(&buffer,sizeof(buffer),1,inm);
    } while (!(buffer == H2T_HEX_COMMA || buffer == H2T_HEX_RIGHTPAR));

    pos2 = ftell(inm);
    
    H2T_SAFE_CALLOC(str,pos2-pos1-1,char,return 0);
    fseek(inm, pos1, 0);
    fread(str, 1, pos2-pos1-1, inm);
    sscanf(str, "%i", &t[dim]);
    H2T_SAFE_FREE(str);
    fread(&buffer,sizeof(buffer),1,inm);
    dim += 1;
  }

  if( !(dim == 3) ) {
    fprintf(stderr,"  ** ERROR : WRONG DIMENSION FOR INPUT DATA."
                   " EXPECTED DIMENSION = 3. PROVIDED DIMENSION = %i \n",dim);
    return -1;
  }

  /* Read data type : look for key 'descr' in header dictionnary */
  fseek(inm,0,0);
  H2T_SAFE_CALLOC(str,5,char,return 0);

  while (dataSize == 0) {
    while (!(buffer == H2T_HEX_APOS)) {
      fread(&buffer,sizeof(buffer),1,inm);
    }

    fread(&str[i],sizeof(char),5,inm);
    if (!strcmp(str,"descr")) {
      /* read header until type-specifying integers are met */
      while (!((buffer >= H2T_HEX_ZERO) && (buffer <= H2T_HEX_NINE))) {
        fread(&buffer,sizeof(buffer),1,inm);
      }
      /* read data size */
      sscanf(&buffer, "%lu", &dataSize);
      H2T_SAFE_FREE(str);
    }
  }

  /* Reach end of header */
  while (!(buffer == H2T_HEX_LINE_FEED)) {
    fread(&buffer,sizeof buffer,1,inm);
  }

  np   = (t[0]+1)*(t[1]+1)*(t[2]+1);
  nhex = t[0]*t[1]*t[2];
  ne   = 4 * (t[0]+t[1]+t[2]);

  if ( H2T_Set_meshSize(mmgMesh,np,nhex,0,ne) != 1 ) {
    return -1;
  }

  H2T_SAFE_CALLOC(*tabhex,9*(nhex+1),int,return 0);

  /* Point coordinates from grid indices */
  ref = 0;
  pos = 0;
  fprintf(stdout,"  READING %" MMG5_PRId " VERTICES\n",np);
  for (i=0;i<=t[0];i++) {
    for (j=0;j<=t[1];j++) {
      for (k=0;k<=t[2];k++) {
        if ( H2T_Set_vertex(mmgMesh,(double)i ,(double)j ,(double)k,ref,++pos) != 1 ) {
          return -1;
        }
      } 
    }
  }

  /* Hexahedra */
  pos = 1;
  fprintf(stdout,"  READING %" MMG5_PRId " HEXA\n",nhex);
  for ( i=0; i<t[0]; ++i ) {
    for ( j=0; j<t[1]; ++j ) {
      for ( k=0; k<t[2]; ++k ) {
        MMG5_int iadr = 9*pos;

        /* Hexa vertices */
        (*tabhex)[iadr+0] = H2T_npy_point_index(i  ,j  ,k  ,t);
        (*tabhex)[iadr+1] = H2T_npy_point_index(i+1,j  ,k  ,t);
        (*tabhex)[iadr+2] = H2T_npy_point_index(i+1,j+1,k  ,t);
        (*tabhex)[iadr+3] = H2T_npy_point_index(i  ,j+1,k  ,t);
        (*tabhex)[iadr+4] = H2T_npy_point_index(i  ,j  ,k+1,t);
        (*tabhex)[iadr+5] = H2T_npy_point_index(i+1,j  ,k+1,t);
        (*tabhex)[iadr+6] = H2T_npy_point_index(i+1,j+1,k+1,t);
        (*tabhex)[iadr+7] = H2T_npy_point_index(i  ,j+1,k+1,t);

        /* Hexa references */
        fread(&(*tabhex)[iadr+8],dataSize,1,inm);
        ++pos;
      }
    }
  }

  /* Edges */
  ref = 0;
  pos = 0;
  fprintf(stdout,"  READING %" MMG5_PRId " EDGES\n",ne);
  for (i=0;i<t[0]; ++i) {
    MMG5_int np0, np1;
    np0 = H2T_npy_point_index(i  ,0     ,0  ,t);
    np1 = H2T_npy_point_index(i+1,0     ,0  ,t);
    if (H2T_Set_edge(mmgMesh,np0,np1,ref,++pos) != 1) {
      return -1;
    }

    np0 = H2T_npy_point_index(i  ,t[1],0  ,t);
    np1 = H2T_npy_point_index(i+1,t[1],0  ,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(i  ,0     ,t[2],t);
    np1 = H2T_npy_point_index(i+1,0     ,t[2],t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(i  ,t[1],t[2],t);
    np1 = H2T_npy_point_index(i+1,t[1],t[2],t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);
  }

  for (j=0;j<t[1]; ++j) {
    MMG5_int np0, np1;
    np0 = H2T_npy_point_index(0  ,j  ,0  ,t);
    np1 = H2T_npy_point_index(0  ,j+1,0  ,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(t[0],j  ,0  ,t);
    np1 = H2T_npy_point_index(t[0],j+1,0  ,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(0  ,j  ,t[2],t);
    np1 = H2T_npy_point_index(0  ,j+1,t[2],t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(t[0],j  ,t[2],t);
    np1 = H2T_npy_point_index(t[0],j+1,t[2],t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);
  }

  for (k=0;k<t[2]; ++k) {
    MMG5_int np0, np1;
    np0 = H2T_npy_point_index(0  ,0  ,k  ,t);
    np1 = H2T_npy_point_index(0  ,0  ,k+1,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(t[0],0  ,k  ,t);
    np1 = H2T_npy_point_index(t[0],0  ,k+1,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(0  ,t[1],k  ,t);
    np1 = H2T_npy_point_index(0  ,t[1],k+1,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);

    np0 = H2T_npy_point_index(t[0],t[1],k  ,t);
    np1 = H2T_npy_point_index(t[0],t[1],k+1,t);
    H2T_Set_edge(mmgMesh,np0,np1,ref,++pos);
  }

  /* Ridges */
  fprintf(stdout,"  READING %" MMG5_PRId " RIDGES\n",ne);
  for (i=1;i<=ne;i++) {
    MMG3D_Set_ridge(mmgMesh,i);
  }

  fclose(inm);

  mmgMesh->ne = 0;
  mmgMesh->nenil = 0;
  for (k=mmgMesh->nenil; k<mmgMesh->nemax-1; k++)
    mmgMesh->tetra[k].v[3] = k+1;

  return nhex;
}

int H2T_loadMesh(MMG5_pMesh mmgMesh,int** tabhex,char *filename) {
  FILE*            inm;
  char             data[128],chaine[128];
  double           x,y,z;
  int              dim,np,nhex,ref,k,iadr,na,nr;
  long             posnp,posna,posnhex,posnr;
  MMG5_int         v0, v1;

  posnp = posna = posnhex = posnr = 0;
  np = nhex = na = nr = 0;

  strcpy(data,filename);
  if( !(inm = fopen(data,"r")) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",data);
    return 0;
  }
  fprintf(stdout,"  %%%% %s OPENED\n",data);

  strcpy(chaine,"D");

  while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
    if(!strncmp(chaine,"Dimension",strlen("Dimension"))) {
      fscanf(inm,"%d",&dim);
      if(dim!=3) {
        fprintf(stdout,"BAD DIMENSION : %d\n",dim);
        return -1;
      }
      continue;
    } else if(!strncmp(chaine,"Vertices",strlen("Vertices"))) {
      fscanf(inm,"%d",&np);
      posnp = ftell(inm);
      continue;
    } else if(!strncmp(chaine,"Hexahedra",strlen("Hexahedra"))) {
      fscanf(inm,"%d",&nhex);
      posnhex = ftell(inm);
      continue;
    } else if(!strncmp(chaine,"Edges",strlen("Edges"))) {
      fscanf(inm,"%d",&na);
      posna = ftell(inm);
      continue;
    } else if(!strncmp(chaine,"Ridges",strlen("Ridges"))) {
      fscanf(inm,"%d",&nr);
      posnr = ftell(inm);
      continue;
    }
  }

  if ( H2T_Set_meshSize(mmgMesh,np,nhex,0,na) != 1 ) {
    return -1;
  }
  H2T_SAFE_CALLOC(*tabhex,9*(nhex+1),int,return 0);

  rewind(inm);
  fseek(inm,posnp,SEEK_SET);
  fprintf(stdout,"  READING %d VERTICES\n",np);
  for (k=1; k<=np; k++) {
    fscanf(inm,"%lf %lf %lf %d",&x,&y,&z,&ref);
    if ( H2T_Set_vertex(mmgMesh,x  ,y  ,z  ,ref,  k) != 1 ) {
      return -1;
    }
  }

  fseek(inm,posnhex,SEEK_SET);
  fprintf(stdout,"  READING %d HEXA\n",nhex);
  for (k=1; k<=nhex; k++) {
    iadr = 9*k;
    fscanf(inm,"%d %d %d %d %d %d %d %d %d",&(*tabhex)[iadr+0],&(*tabhex)[iadr+1]
	   ,&(*tabhex)[iadr+2],&(*tabhex)[iadr+3],&(*tabhex)[iadr+4],
	   &(*tabhex)[iadr+5],&(*tabhex)[iadr+6],&(*tabhex)[iadr+7],&(*tabhex)[iadr+8]);
  }

  if(na) {
    fseek(inm,posna,SEEK_SET);
    fprintf(stdout,"  READING %d EDGES\n",na);
    for (k=1; k<=na; k++) {
      fscanf(inm,"%" MMG5_PRId " %" MMG5_PRId " %d", &v0, &v1, &ref);
      if ( H2T_Set_edge(mmgMesh, v0, v1, ref, k) != 1 ) {
	return -1;
      }
    }
  }

  if(nr) {
    fseek(inm,posnr,SEEK_SET);
    fprintf(stdout,"  READING %d RIDGES\n",na);
    for (k=1; k<=nr; k++) {
      fscanf(inm,"%" MMG5_PRId, &v0);
      if ( !MMG3D_Set_ridge(mmgMesh,v0) ) {
	return -1;
      }
    }
  }

  fclose(inm);

  mmgMesh->ne = 0;
  mmgMesh->nenil = 0;
  for (k=mmgMesh->nenil; k<mmgMesh->nemax-1; k++)
    mmgMesh->tetra[k].v[3] = k+1;

  return nhex;
}
