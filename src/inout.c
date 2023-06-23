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
#include <string.h>
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>


static inline
int H2T_npy_point_index(int i,int j,int k, npy_intp *n) {
  return n[1]*n[2]*i + n[2]*j + k + 1;
}

#warning clean findPython and python option in CMAKE
#warning add USE_PYTHON

int H2T_loadNpyArray(MMG5_pMesh mmgMesh,int** tabhex,char *filename) {
  PyObject *pName, *pModule, *pDict, *pFunc,*pValue;
  int np,i,j,k,nbhex;

  Py_Initialize();
  pName = PyUnicode_FromString("numpy");

  /* Error checking of pName left out */
  pModule = PyImport_Import(pName);
  Py_DECREF(pName);

  assert ( pModule );

  pFunc = PyObject_GetAttrString(pModule, "load");
  assert ( pFunc && PyCallable_Check(pFunc) );

  pName = PyUnicode_FromString(filename);
  pValue = PyObject_CallObject(pFunc, pName);

  Py_DECREF(pName);
  Py_DECREF(pFunc);


  if ( !pValue ) {
    /* File reading has failed */
    fprintf(stderr,"  %%%% %s UNABLE TO LOAD FILE\n",filename);
    return -1;
  }

  int dim = PyArray_NDIM(pValue);

  if ( dim != 3 ) {
    fprintf(stderr," ## Error: Unexpected array dimension (%d while instead of"
            " a 3D array)\n",dim);
    Py_DECREF(pValue);
    return -1;
  }

  npy_intp *n = PyArray_DIMS(pValue);

  /* Mesh size evaluation and mesh allocation */
  np = n[0]*n[1]*n[2];
  nbhex = (n[0]-1)*(n[1]-1)*(n[2]-1);

  if ( H2T_Set_meshSize(mmgMesh,np,nbhex,0,0) != 1 ) {
    Py_DECREF(pValue);
    return -1;
  }
  *tabhex = (int*) malloc(9*(nbhex+1)*sizeof(int));

  /* Point coordinates from grid indices */
  int pos = 0;
  int ref = 0;
  for ( i=0; i<n[0]; ++i ) {
    for ( j=0; j<n[1]; ++j ) {
      for ( k=0; k<n[2]; ++k ) {
        if ( H2T_Set_vertex(mmgMesh,(double)i ,(double)j ,(double)k,ref,++pos) != 1 ) {
          return -1;
        }
      }
    }
  }

  /* Hexahedra */
  pos = 1;
  for ( i=0; i<n[0]; ++i ) {
    for ( j=0; j<n[1]; ++j ) {
      for ( k=0; k<n[2]; ++k ) {
        int iadr = 9*pos;

        /* Hexa vertices */
        *tabhex[iadr+0] = H2T_npy_point_index(i  ,j  ,k  ,n);
        *tabhex[iadr+1] = H2T_npy_point_index(i+1,j  ,k  ,n);
        *tabhex[iadr+2] = H2T_npy_point_index(i+1,j+1,k  ,n);
        *tabhex[iadr+3] = H2T_npy_point_index(i  ,j+1,k  ,n);
        *tabhex[iadr+4] = H2T_npy_point_index(i  ,j  ,k+1,n);
        *tabhex[iadr+5] = H2T_npy_point_index(i+1,j  ,k+1,n);
        *tabhex[iadr+6] = H2T_npy_point_index(i+1,j+1,k+1,n);
        *tabhex[iadr+7] = H2T_npy_point_index(i  ,j+1,k+1,n);

        /* Hexa references */
        int *tmp;
        tmp = (int*)PyArray_GETPTR3(pValue,i,j,k);

        *tabhex[iadr+7] = *tmp;
        ++pos;
      }
    }
  }
  Py_DECREF(pValue);

  mmgMesh->ne = 0;
  mmgMesh->nenil = 0;
  for (k=mmgMesh->nenil; k<mmgMesh->nemax-1; k++)
    mmgMesh->tetra[k].v[3] = k+1;

  return nbhex;

}


int H2T_loadMesh(MMG5_pMesh mmgMesh,int* tabhex,int nbhex,char *filename) {
  FILE*            inm;
  char             data[128],chaine[128];
  double           x,y,z;
  int              dim,np,nhex,ref,k,iadr,na,nr;
  long             posnp,posna,posnhex,posnr;
  MMG5_int         v0, v1;

  posnp = posna = posnhex = posnr = 0;

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

  if ( H2T_Set_meshSize(mmgMesh,np,nbhex,0,na) != 1 ) {
    return -1;
  }

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
    fscanf(inm,"%d %d %d %d %d %d %d %d %d",&tabhex[iadr+0],&tabhex[iadr+1]
	   ,&tabhex[iadr+2],&tabhex[iadr+3],&tabhex[iadr+4],
	   &tabhex[iadr+5],&tabhex[iadr+6],&tabhex[iadr+7],&tabhex[iadr+8]);
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
