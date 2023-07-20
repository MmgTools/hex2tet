/**
 * \file libhex2tet.h
 * \brief C API for hex2tet library.
 * \author C. Dobrzynski (Bx INP/Inria/UBordeaux)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 *
 */

#ifndef _H2TLIB_H
#define _H2TLIB_H

#ifdef __cplusplus
extern "C" {
#endif

#include "mmg/mmg3d/libmmg3d.h"

#define  MG_NUL       (1 << 14)

/**
 * \def H2T_SUCCESS
 *
 * Return value for success.
 *
 */
#define H2T_SUCCESS       0
/**
 * \def H2T_LOWFAILURE
 *
 * Return value if the remesh process failed but we can save a conform
 * mesh.
 *
 */
#define H2T_LOWFAILURE    1
/**
 * \def H2T_STRONGFAILURE
 *
 * Return value if the remesh process failed and the mesh is
 * non-conform.
 *
 */
#define H2T_STRONGFAILURE 2
/**
 * \def H2T_ARG_hexa
 *
 * Pointer towards hexahedra array.
 *
 * \remark we cannot use an enum because used in
 * variadic functions).
 */
#define H2T_ARG_phexa  11
/**
 * \param mmgMesh mesh structure with only vertices.
 * \param hexa tab of hexahedra (size 9*(nbhexa+1) : 8 vertices per hexa + one ref).
 * \param nbhexa number of hexa
 *
 * \return \ref H2T_SUCCESS if success.
 * \return \ref H2T_LOWFAILURE if failed but a conform mesh is saved.
 * \return \ref H2T_STRONGFAILURE if failed and we can't save the mesh.
 *
 * main library function: cut hexa into tetra.
 *
 * \remark Fortran interface:
 * >   SUBROUTINE H2T_LIBHEX2TET(mmgMesh,hexa,nbHexa,filename,strlen0,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)     :: mmgMesh\n
 * >     INTEGER, DIMENSION(*), INTENT(IN) :: hexa\n
 * >     INTEGER, INTENT(IN)               :: nbHexa\n
 * >     CHARACTER(LEN=*), INTENT(IN)      :: filename\n
 * >     INTEGER, INTENT(IN)               :: strlen0\n
 * >     INTEGER, INTENT(OUT)              :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int H2T_libhex2tet(MMG5_pMesh mmgMesh,int* hexa,MMG5_int nbhexa );

/**
 * \param mmgMesh pointer toward the mesh.
 * \param tabhex array of hexa.
 * \param filname name of input file.
 *
 * \return  0 if the file is not found, -1 if we detect mismatch parameters or we
 * fail, 1 otherwise.
 *
 * Read mesh data.
 *
 * * \remark Fortran interface:
 * >   SUBROUTINE H2T_LOADMESH(mmgMesh,tabhex,nbhex,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)     :: mmgMesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*)  :: tabhex\n
 * >     INTEGER, INTENT(OUT)              :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int H2T_loadMesh(MMG5_pMesh mmgMesh,int** tabhex,char *filename);

/**
 * \param mmgMesh pointer toward the mesh.
 * \param tabhex pointer towards array of hexa.
 * \param filename name of input file.
 *
 * \return  0 if the file is not found, -1 if we detect mismatch parameters or we
 * fail. If successful, number of read hexahedra.
 *
 * Read mesh data from .npy file
 *
 * * \remark Fortran interface:
 * >   SUBROUTINE H2T_LOADNPY(mmgMesh,tabhex,filename,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT)     :: mmgMesh\n
 * >     INTEGER(MMG5F_INT), DIMENSION(*)  :: tabhex\n
 * >     CHARACTER(LEN=*), INTENT(IN)      :: filename\n
 * >     INTEGER, INTENT(OUT)              :: retval\n
 * >   END SUBROUTINE\n
 *
 */
int H2T_loadNpy(MMG5_pMesh mmgMesh,int** tabhex,char *filename);

/**
 * \param starter dummy argument used to initialize the variadic argument
 * list
 * \param ... variadic arguments that depend to the library function that you
 * want to call.
 *
 * For now, you need to provide the following arguments :
 * MMG5_ARG_start, to speficy the beginning of the variadic arg list
 * MMG5_ARG_ppMesh, to say that you will provide your mesh structure
 * &your_mesh, pointer toward your mesh structure
 * MMG5_ARG_ppMet or MMG5_ARG_ppLs or MMG5_ARG_ppDisp to say that you will give
 * a MMG5_pSol structure to strore (resp.) a metric, a level-set or a displacement
 * &your_sol_stuct, a pointer toward your solution structure
 * MMG5_ARG_end  , to speficy the end of the variadic arg list
 *
 * example:
 * H2T_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh, &your_mesh, MMG5_ARG_ppMet
 * MMG5_ARG_ppMet, &your_metric,MMG5_ARG_end).
 *
 * \return 1 if success, 0 if fail
 *
 * Wrapper for the MMG3D_Init_mesh function (MMG structures allocation and initialization).
 *
 * \remark No fortran interface to allow variadic arguments.
 *
 */
  int H2T_Init_mesh(const int starter,...);

/**
 * \param mesh pointer toward the mesh structure.
 * \param np number of vertices.
 * \param nhexa number of hexahedra.
 * \return 0 if failed, 1 otherwise.
 *
 * Set the number of vertices and of hexahedra of the mesh
 *
 * \remark Fortran interface:
 * >   SUBROUTINE H2T_SET_MESHSIZE(mesh,np,nhexa,nquad,na,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     INTEGER                       :: np,nhexa,nquad,na\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  H2T_Set_meshSize(MMG5_pMesh mesh,int np,int nhexa,int nquad,int na);

/**
 * \param mesh pointer toward the mesh structure.
 * \param c0 coordinate of the point along the first dimension.
 * \param c1 coordinate of the point along the second dimension.
 * \param c2 coordinate of the point along the third dimension.
 * \param ref point reference.
 * \param pos position of the point in the mesh.
 * \return 1.
 *
 * Set vertex of coordinates \a c0, \a c1,\a c2 and reference \a ref
 * at position \a pos in mesh structure
 *
 * \remark Fortran interface:
 * >   SUBROUTINE H2T_SET_VERTEX(mesh,c0,c1,c2,ref,pos,retval)\n
 * >     MMG5_DATA_PTR_T,INTENT(INOUT) :: mesh\n
 * >     REAL(KIND=8), INTENT(IN)      :: c0,c1,c2\n
 * >     INTEGER, INTENT(IN)           :: ref,pos\n
 * >     INTEGER, INTENT(OUT)          :: retval\n
 * >   END SUBROUTINE\n
 *
 */
  int  H2T_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                      double c2, int ref,int pos);

  /**
   * \param mesh pointer toward the mesh structure.
   * \param v0 first end point of the edge.
   * \param v1 second end point of the edge.
   * \param ref edge reference.
   * \param pos position of the edge in the mesh.
   * \return 1.
   *
   * Set edge of end points \a v0, \a v1 and reference \a ref
   * at position \a pos in mesh structure
   *
   * \remark Fortran interface:
   * >   SUBROUTINE H2T_SET_EDGE(mesh,v0,v1,ref,pos,retval)\n
   * >     MMG5_DATA_PTR_T,INTENT(INO UT) :: mesh\n
   * >     INTEGER(MMG5F_INT), INTENT(IN) :: v0,v1\n
   * >     INTEGER, INTENT(IN)            :: ref,pos\n
   * >     INTEGER, INTENT(OUT)           :: retval\n
   * >   END SUBROUTINE\n
   *
   */
  int H2T_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1,
		   MMG5_int ref, MMG5_int pos);

/**
 * \param starter dummy argument used to initialize the variadic argument
 * list.
 * \param ... variadic arguments that depend to the library function that you
 * have call.
 *
 * \return 1 if success, 0 if fail
 *
 * Deallocations before return.
 *
 * \remark we pass the structures by reference in order to have argument
 * compatibility between the library call from a Fortran code and a C code.
 *
 * \remark no Fortran interface to allow variadic args.
 *
 */
int H2T_Free_all(const int starter,...);


#ifdef __cplusplus
}
#endif

#endif
