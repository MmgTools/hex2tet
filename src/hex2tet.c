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
#include <string.h>
#include <stdlib.h>
#include <signal.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <limits.h>

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
 * \param fmt file format.
 *
 * \return The name of the file format in a string.
 *
 * Print the name of the file format associated to \a fmt.
 *
 */
static inline
const char* H2T_Get_formatName(enum MMG5_Format fmt)
{
  switch (fmt)
  {
  case H2T_FMT_MeditASCII:
    return "H2T_FMT_MeditASCII";
    break;
  case H2T_FMT_Numpy:
    return "H2T_FMT_Numpy";
    break;
  default:
    return "H2T_Unknown";
  }
}

/**
 * \param filename string containing a filename
 *
 * \return pointer toward the filename extension or toward the end of the string
 * if no extension have been founded
 *
 * Get the extension of the filename string. Do not consider '.o' as an extension.
 *
 * \warning copy past from Mmg
 */
static inline
char *H2T_Get_filenameExt( char *filename ) {
  const char pathsep='/';
  char       *dot,*lastpath;

  if ( !filename ) {
    return NULL;
  }

  dot = strrchr(filename, '.');
  lastpath = (pathsep == 0) ? NULL : strrchr (filename, pathsep);

  if ( (!dot) || dot == filename || (lastpath>dot) || (!strcmp(dot,".o")) ) {
    /* No extension */
    return filename + strlen(filename);
  }

  return dot;
}

/**
 * \param ptr pointer toward the file extension (dot included)
 * \param fmt default file format.
 *
 * \return and index associated to the file format detected from the extension.
 *
 * Get the wanted file format from the mesh extension. If \a fmt is provided, it
 * is used as default file format (\a ptr==NULL), otherwise, the default file
 * format is the medit one.
 *
 */
int H2T_Get_format( char *ptr, int fmt ) {
  /* Default is the Medit file format or a format given as input */
  int defFmt = fmt;

  if ( !ptr ) return defFmt;

  if ( !strncmp( ptr,".mesh",strlen(".mesh") ) ) {
    return H2T_FMT_MeditASCII;
  }
  else if ( !strncmp ( ptr,".npy",strlen(".npy") ) ) {
    return H2T_FMT_Numpy;
  }

  return defFmt;
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
  char            chaine[128],*ptr;
  int             *hexa,nbhex,ier,fmtin;

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

  ptr   = H2T_Get_filenameExt(mmgMesh->namein);
  fmtin = H2T_Get_format(ptr,MMG5_FMT_MeditASCII);

  switch ( fmtin ) {

    case ( H2T_FMT_Numpy ):

      nbhex = H2T_loadNpy(mmgMesh,&hexa,mmgMesh->namein);
      break;

    case ( H2T_FMT_MeditASCII ):

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
          H2T_SAFE_CALLOC(hexa,9*(nbhex+1),int,return 0);
          assert(hexa);
          break;
        }
      }
      fclose(inm);

      nbhex = H2T_loadMesh(mmgMesh,hexa,nbhex,mmgMesh->namein);

      break;

    default:
      fprintf(stderr,"  ** I/O AT FORMAT %s NOT IMPLEMENTED.\n",H2T_Get_formatName(fmtin) );
      return H2T_STRONGFAILURE;
  }

  if ( nbhex < 0 ) {
    return H2T_STRONGFAILURE;
  }

  /** call hex2tet library */
  ier = H2T_libhex2tet(mmgMesh,hexa,nbhex);

  MMG3D_saveMesh(mmgMesh,mmgMesh->nameout);

  /** free structures */
  H2T_SAFE_FREE(hexa);

  MMG3D_Free_all(MMG5_ARG_start,
                 MMG5_ARG_ppMesh,&mmgMesh,MMG5_ARG_ppMet,&mmgSol,
                 MMG5_ARG_end);

  return ier;
}
