/**
 * \file API_functionsf_h2t.c
 * \brief Fortran API functions definitions for the hex2tet library.
 *
 * \author Algiane Froehly (InriaSoft)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

#include <stdlib.h>

/**
 * See \ref H2T_Init_mesh function in libhex2tet.h file.
 */
FORTRAN_VARIADIC ( H2T_INIT_MESH, h2t_init_mesh,
                   (const int starter, ... ),
                   va_list argptr;
                   int     ier;

                   va_start(argptr, starter);

                   ier = _H2T_Init_mesh_var(argptr);

                   va_end(argptr);

                   if ( !ier ) exit(EXIT_FAILURE);

                   return;
  )

/**
 * See \ref H2T_Free_all function in libhex2tet.h file.
 */
FORTRAN_VARIADIC ( H2T_FREE_ALL, h2t_free_all,
                   (const int starter, ... ),
                   va_list argptr;
                   int     ier;

                   va_start(argptr, starter);

                   ier = H2T_Free_all_var(argptr);

                   va_end(argptr);

                   if ( !ier ) exit(EXIT_FAILURE);

                   return;
  )

/**
 * See \ref H2T_Set_meshSize function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_SET_MESHSIZE,h2t_set_meshsize,
             (MMG5_pMesh *mesh, int *np, int *nhexa,
              int *nquad, int *na, int *retval),
             (mesh,np,nhexa,nquad,na,retval)) {
  *retval = H2T_Set_meshSize(*mesh,*np,*nhexa,*nquad,*na);
  return;
}

/**
 * See \ref H2T_loadMesh function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_LOADMESH,h2t_loadmesh,
             (MMG5_pMesh *mesh, int **tabhex,
              char* filename, int *strlen0,int *retval),
             (mesh,tabhex,filename,strlen0,retval)) {

  *retval = H2T_loadMesh(*mesh,tabhex,filename);
  return;
}

/**
 * See \ref H2T_loadNpy function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_LOADNPY,h2t_loadnpy,
             (MMG5_pMesh *mesh, int **tabhex,
              char* filename, int *retval),
             (mesh,tabhex,filename, retval)) {

  *retval = H2T_loadNpy(*mesh,tabhex,filename);
  return;
}

/**
 * See \ref H2T_Set_vertex function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_SET_VERTEX,h2t_set_vertex,
             (MMG5_pMesh *mesh, double* c0, double* c1, double* c2, int* ref,
              int* pos, int* retval),
             (mesh,c0,c1,c2,ref,pos,retval)) {

  *retval = H2T_Set_vertex(*mesh,*c0,*c1,*c2,*ref,*pos);
  return;
}

/**
 * See \ref H2T_Set_edge function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_SET_EDGE,h2t_set_edge,
             (MMG5_pMesh *mesh, MMG5_int* v0, MMG5_int* v1, MMG5_int* ref,
              MMG5_int* pos, int* retval),
             (mesh,v0,v1,ref,pos,retval)) {

  *retval = H2T_Set_edge(*mesh,*v0,*v1,*ref,*pos);
  return;
}

/**
 * See \ref H2T_Set_hexahedron function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_SET_HEXAHEDRON,h2t_set_hexahedron,
             (int *hexTab,
              int* i0,int* i1,int* i2,int* i3,int* i4,int* i5,int* i6,int* i7,
              int* ref,int* pos,int* retval),
             (hexTab,i0,i1,i2,i3,i4,i5,i6,i7,ref,pos,retval)) {

  *retval = H2T_Set_hexahedron(hexTab,*i0,*i1,*i2,*i3,*i4,*i5,*i6,*i7,*ref,*pos);
  return;
}

/**
 * See \ref H2T_Set_quadrilateral function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_SET_QUADRILATERAL,h2t_set_quadrilateral,
             (MMG5_pMesh *mesh,
              int* i0,int* i1,int* i2,int* i3,
              int* ref,int* pos,int* retval),
             (mesh,i0,i1,i2,i3,ref,pos,retval)) {

  *retval = H2T_Set_quadrilateral(*mesh,*i0,*i1,*i2,*i3,*ref,*pos);
  return;
}

/**
 * See \ref H2T_libhex2tet function in \ref libhex2tet.h file.
 */
FORTRAN_NAME(H2T_LIBHEX2TET,h2t_libhex2tet,
             (MMG5_pMesh *mmgMesh, int **hexa, int* nbHexa, int* retval),
             (mmgMesh,hexa,nbHexa,retval)) {

  *retval = H2T_libhex2tet(*mmgMesh,hexa,*nbHexa);
  return;
}
