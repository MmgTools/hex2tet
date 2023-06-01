/**
 * \file API_functions_h2t.c
 * \brief C API functions definitions for the hex2tet library.
 *
 * \author Cecile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>


int H2T_Init_mesh(const int starter,...) {
  va_list        argptr;
  int            ier;

  va_start(argptr, starter);

  ier = _H2T_Init_mesh_var(argptr);

  va_end(argptr);

  return ier;
}

int _H2T_Init_mesh_var(va_list argptr) {
  MMG5_pMesh     *mesh;
  MMG5_pSol      *sol,*disp;
  int            typArg,meshCount,solCount,metCount,dispCount,ier;

  mesh = NULL;
  disp = sol = NULL;

  meshCount = solCount = dispCount = metCount = 0;

  while ( (typArg = va_arg(argptr,int)) != MMG5_ARG_end )
  {
    switch ( typArg )
    {
    case(MMG5_ARG_ppMesh):
      mesh = va_arg(argptr,MMG5_pMesh*);
      ++meshCount;
      break;
    case(MMG5_ARG_ppMet):
      sol = va_arg(argptr,MMG5_pSol*);
      ++metCount;
      break;
    case MMG5_ARG_ppLs:
      sol = va_arg(argptr,MMG5_pSol*);
      ++solCount;
      break;
    case(MMG5_ARG_ppDisp):
      disp = va_arg(argptr,MMG5_pSol*);
      ++dispCount;
      break;
    default:
      fprintf(stderr,"\n  ## Error: %s:\n"
              " unexpected argument type: %d\n",__func__,typArg);
      fprintf(stderr," Argument type must be one"
              " of the MMG5_ARG* preprocessor variable:"
              " MMG5_ARG_ppMesh, MMG5_ARG_ppMet,"
              "  MMG5_ARG_ppLs, MMG5_ARG_ppDisp\n");
      return 0;
    }
  }

  if ( meshCount !=1 ) {
    fprintf(stderr,"\n  ## Error: %s:\n"
            " you need to initialize the mesh structure that"
            " will contain your mesh.\n",__func__);
    return 0;
  }
  if ( !sol ) {
    fprintf(stderr,"\n  ## Error: %s:\n"
            " you need to initialize a solution structure"
            " (of type MMG5_pSol and indentified by the MMG5_ARG_ppMet or"
            " MMG5_ARG_ppLs preprocessor variable) that will contain the output"
            " mesh metric"
            " informations, and the input one, if provided.\n.",__func__);
    return 0;
  }


  if ( solCount + metCount > 1 ) {
    fprintf(stderr,"\n  ## Error: %s:\n"
            " you can provide only one argument of type %s, %s or %s\n",
            __func__,"MMG5_ARG_ppMet","MMG5_ARG_ppLs","MMG5_ARG_ppDisp");
    return 0;
  }

  ier = 0;
  if ( metCount ) {
    assert ( !solCount );
    if ( dispCount ) {
      ier = MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,mesh,
                            MMG5_ARG_ppMet,sol,MMG5_ARG_ppDisp,disp,MMG5_ARG_end);
    }
    else {
      ier = MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,mesh,
                            MMG5_ARG_ppMet,sol,MMG5_ARG_end);
    }
  }
  else if ( solCount ) {
    assert ( !metCount );
    ier = MMG3D_Init_mesh(MMG5_ARG_start,MMG5_ARG_ppMesh,mesh,
                          MMG5_ARG_ppLs,sol,MMG5_ARG_end);

  }

  return ier;
}

 int  H2T_Set_meshSize(MMG5_pMesh mesh,int np,int nhexa,int nquad,int na) {
   int ne,k;

   ne = 6*nhexa;

   if(MMG3D_Set_meshSize(mesh, np, 6*nhexa, 0, 0, 0, na) != 1)
     return 0;

   if ( nquad ) {
     printf("  ## Warning: %s: treatment of input quadrangles not yet implemented:\n"
            "%d quadrangles ignored.\n",__func__,nquad);
   }

   /* Set all tetra as unused */
   mesh->nenil = 1;
   for ( k=mesh->nenil; k<mesh->nemax-1; k++)
     mesh->tetra[k].v[3] = k+1;

   return 1;
 }

  int  H2T_Set_vertex(MMG5_pMesh mesh, double c0, double c1,
                      double c2, int ref,int pos) {

    return MMG3D_Set_vertex(mesh,c0,c1,c2,ref,pos);

  }

int H2T_Set_edge(MMG5_pMesh mesh, MMG5_int v0, MMG5_int v1, int ref, int pos) {

  return MMG3D_Set_edge(mesh, v0, v1, ref, pos);

}
