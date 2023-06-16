/**
 * \file hex2tet.c
 * \brief hex2tet application: convert an hexahedral mesh into a tetrahedral one.
 *
 * \author Cecile Dobrzynski (Bx INP/Inria/UBordeaux)
 * \author Algiane Froehly (InriaSoft)
 *
 * \version 1
 * \copyright GNU Lesser General Public License.
 */

#include "hex2tet.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <signal.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>

/**
 * \param mmgMesh pointer toward the mesh.
 * \param tabhex array of hexa.
 * \param nbhex number of hexa.
 * \param filname name of input file.
 *
 * \return  0 if the file is not found, -1 if we detect mismatch parameters or we
 * fail, 1 otherwise.
 *
 * Read mesh data.
 *
 * \warning todo: add the hexa to the MMG5_mesh structure and read it
 * in the MMG3D_loadmesh function
 *
 */

int H2T_loadMesh(MMG5_pMesh mmgMesh,int* tabhex,int nbhex,char *filename) {
  FILE*            inm;
  char             data[128],chaine[128];
  double           x,y,z;
  int              dim,np,nhex,ref,k,iadr,na,nr;
  long             posnp,posna,posnhex,posnr;
  MMG5_int         v0, v1;

  posnp = posna = posnhex = posnr = 0;
  na = nr = 0;

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



/**
 * \param *prog pointer toward the program name.
 *
 * Print help for hex2tet application
 *
 */
static
int H2T_usage(char *prog) {

  fprintf(stdout,"\nUsage: %s [-v [n]] [opts..] filein [fileout]\n",prog);

  fprintf(stdout,"\n** Generic options :\n");
  fprintf(stdout,"-h        Print this message\n");
  fprintf(stdout,"-v [n]    Tune level of verbosity, [-10..10]\n");
  fprintf(stdout,"-m [n]    Set maximal memory size to n Mbytes\n");
  fprintf(stdout,"-d        Turn on debug mode\n");

  fprintf(stdout,"\n**  File specifications\n");
  fprintf(stdout,"-in  file  input triangulation\n");
  fprintf(stdout,"-out file  output triangulation\n");

  return 1;
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \param mesh pointer toward the mesh structure.
 * \param met pointer toward the sol structure.
 * \return 1 if we want to run Mmg after, 0 if not or if fail.
 *
 * Store command line arguments.
 *
 */
static
int H2T_parsar(int argc,char *argv[],MMG5_pMesh mesh,MMG5_pSol met) {
  int     i;
  char    namein[128];

  /* read arguments. */
  i = 1;
  while ( i < argc ) {
    if ( *argv[i] == '-' ) {
      switch(argv[i][1]) {
      case '?':
        H2T_usage(argv[0]);
        return 0;
      case 'd':
        /* debug */
        if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_debug,1) ) {
          return 0;
        }
        break;
      case 'h':
        H2T_usage(argv[0]);
        return 0;
        break;
      case 'i':
        if ( !strcmp(argv[i],"-in") ) {
          if ( ++i < argc && isascii(argv[i][0]) && argv[i][0]!='-') {
            if ( !MMG3D_Set_inputMeshName(mesh, argv[i]) )
              return 0;

          }else{
            fprintf(stderr,"Missing filname for %c%c\n",argv[i-1][1],argv[i-1][2]);
            H2T_usage(argv[0]);
            return 0;
          }
        }
        break;
      case 'm':  /* memory */
        if ( ++i < argc && isdigit(argv[i][0]) ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_mem,atoi(argv[i])) )
            return 0;
        }
        else {
          fprintf(stderr,"Missing argument option %c\n",argv[i-1][1]);
          H2T_usage(argv[0]);
          return 0;
        }
        break;
      case 'o':
        if ( ++i < argc && isascii(argv[i][0])  && argv[i][0]!='-') {
          if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
            return 0;
        }else{
          fprintf(stderr,"Missing filname for %c%c%c\n",
                  argv[i-1][1],argv[i-1][2],argv[i-1][3]);
          H2T_usage(argv[0]);
          return 0;
        }
        break;
      case 'v':
        if ( ++i < argc ) {
          if ( isdigit(argv[i][0]) ||
               (argv[i][0]=='-' && isdigit(argv[i][1])) ) {
            if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,atoi(argv[i])) )
              return 0;
          }
          else {
            i--;
            fprintf(stderr,"Missing argument option %s\n",argv[i]);
          }
        }
        else {
          fprintf(stderr,"Missing argument option %s\n",argv[i-1]);
          H2T_usage(argv[0]);
          return 0;
        }
        break;
      default:
        fprintf(stderr,"Unrecognized option %s\n",argv[i]);
        H2T_usage(argv[0]);
        return 0;
      }
    }
    else {
      if ( mesh->namein == NULL ) {
        if ( !MMG3D_Set_inputMeshName(mesh,argv[i]) )
          return 0;
        if ( mesh->info.imprim == -99 ) {
          if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,5) )
            return 0;
        }
      }
      else if ( mesh->nameout == NULL ) {
        if ( !MMG3D_Set_outputMeshName(mesh,argv[i]) )
          return 0;
      }
      else {
        fprintf(stdout,"Argument %s ignored\n",argv[i]);
        H2T_usage(argv[0]);
        return 0;
      }
    }
    i++;
  }

  /* check file names */
  if ( mesh->info.imprim == -99 ) {
    fprintf(stdout,"\n  -- PRINT (0 10(advised) -10) ?\n");
    fflush(stdin);
    fscanf(stdin,"%d",&i);
    if ( !MMG3D_Set_iparameter(mesh,met,MMG3D_IPARAM_verbose,i) )
      return 0;
  }

  if ( mesh->namein == NULL ) {
    fprintf(stdout,"  -- INPUT MESH NAME ?\n");
    fflush(stdin);
    fscanf(stdin,"%127s",namein);
    if ( !MMG3D_Set_inputMeshName(mesh,namein) )
      return 0;
  }

  if ( mesh->nameout == NULL ) {
    if ( !MMG3D_Set_outputMeshName(mesh,"") )
      return 0;
  }

  if ( met->namein == NULL ) {
    if ( !MMG3D_Set_inputSolName(mesh,met,"") )
      return 0;
  }
  if ( met->nameout == NULL ) {
    if ( !MMG3D_Set_outputSolName(mesh,met,"") )
      return 0;
  }

  return 1;
}

/**
 * \param argc number of command line arguments.
 * \param argv command line arguments.
 * \return \ref H2T_SUCCESS if success.
 * \return \ref H2T_LOWFAILURE if failed but a conform mesh is saved.
 * \return \ref H2T_STRONGFAILURE if failed and we can't save the mesh.
 *
 * Main program to convert an hexa mesh into a tetra one.
 *
 */
int main(int argc,char *argv[]) {
  FILE*           inm;
  MMG5_pMesh      mmgMesh;
  MMG5_pSol       mmgSol;
  char            chaine[128];
  int             *hexa,nbhex,ier;

  fprintf(stdout,"\n  -- H2T, Release %s (%s) \n",H2T_VER,H2T_REL);
  fprintf(stdout,"     %s\n",H2T_CPY);
  fprintf(stdout,"     %s %s\n\n",__DATE__,__TIME__);

  /** Mesh initialization */
  mmgMesh = NULL;
  mmgSol  = NULL;
  hexa    = NULL;

  MMG3D_Init_mesh(MMG5_ARG_start,
                  MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                  MMG5_ARG_end);

  /* reset default values for file names */
  if ( !MMG3D_Free_names(MMG5_ARG_start,
                         MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                         MMG5_ARG_end) )
    return H2T_STRONGFAILURE;

  /* command line */
  if ( !H2T_parsar(argc,argv,mmgMesh,mmgSol) )  return H2T_STRONGFAILURE;

  /* Input data and creation of the hexa array */
  if( !(inm = fopen(mmgMesh->namein,"r")) ) {
    fprintf(stderr,"  ** %s  NOT FOUND.\n",mmgMesh->namein);
    return 0;
  }

  nbhex = 0;
  strcpy(chaine,"D");
  while(fscanf(inm,"%s",&chaine[0])!=EOF && strncmp(chaine,"End",strlen("End")) ) {
    if(!strncmp(chaine,"Hexahedra",strlen("Hexahedra"))) {
      fscanf(inm,"%d",&nbhex);
      fprintf(stdout,"  READING %d HEXA\n",nbhex);
      hexa = (int*) malloc(9*(nbhex+1)*sizeof(int));
      assert(hexa);
      break;
    }
  }
  fclose(inm);

  nbhex = H2T_loadMesh(mmgMesh,hexa,nbhex,mmgMesh->namein);

  if ( nbhex < 0 ) {
    return H2T_STRONGFAILURE;
  }

  /** call hex2tet library */
  ier = H2T_libhex2tet(mmgMesh,hexa,nbhex);

  MMG3D_saveMesh(mmgMesh,mmgMesh->nameout);

  /** free structures */
  free(hexa);

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return ier;
}
