CMAKE_MINIMUM_REQUIRED(VERSION 3.5.0)

PROJECT (hex2tet)

INCLUDE(CMakeDependentOption)
INCLUDE(cmake/modules/macros.cmake)

# Executable path
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
SET(LIBRARY_OUTPUT_PATH    ${PROJECT_BINARY_DIR}/lib)


SET(H2T_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
SET(H2T_BUILD_DIR  ${PROJECT_BINARY_DIR}/src)
SET(H2T_INCLUDE    ${PROJECT_BINARY_DIR}/include/${PROJECT_NAME} )


FILE ( MAKE_DIRECTORY  ${H2T_BUILD_DIR} )


###############################################################################
#####
#####         Release version and date
#####
###############################################################################
SET (CMAKE_RELEASE_VERSION_MAJOR "1")
SET (CMAKE_RELEASE_VERSION_MINOR "1")
SET (CMAKE_RELEASE_VERSION_PATCH "0")
SET (CMAKE_RELEASE_DATE "Nov 5, 2024")

SET (CMAKE_RELEASE_VERSION
  "${CMAKE_RELEASE_VERSION_MAJOR}.${CMAKE_RELEASE_VERSION_MINOR}.${CMAKE_RELEASE_VERSION_PATCH}")

# Create hex2tet.h file with the good release infos.
CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/src/hex2tet.h.in
  ${H2T_BUILD_DIR}/hex2tet.h @ONLY)
INCLUDE_DIRECTORIES ( ${H2T_BUILD_DIR} )

###############################################################################
#####
#####         Main CMake compilation variables
#####
###############################################################################

# To see flags and options of compilation
#SET(CMAKE_VERBOSE_MAKEFILE TRUE)

# Find "exotic" compilers
IF (${CMAKE_C_COMPILER} MATCHES Clang OR ${CMAKE_C_COMPILER} MATCHES clang)
  # using clang
  SET(CMAKE_COMPILER_IS_CLANG TRUE)
ELSEIF(${CMAKE_C_COMPILER} MATCHES Icc OR ${CMAKE_C_COMPILER} MATCHES icc)
  # using icc
  SET(CMAKE_COMPILER_IS_INTEL TRUE)
ENDIF()

# Hide some options
MARK_AS_ADVANCED(CMAKE_OSX_ARCHITECTURES
  CMAKE_OSX_DEPLOYMENT_TARGET
  CMAKE_OSX_SYSROOT)

IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS " -Wno-char-subscripts ${CMAKE_C_FLAGS}")
  IF(APPLE)
    # Add flags to the compiler to work on old mac
    ADD_DEFINITIONS( -mmacosx-version-min=10.4 -arch x86_64)

    # To avoid pbs with binary files...
    SET(CMAKE_EXE_LINKER_FLAGS "-arch x86_64 ${CMAKE_EXE_LINKER_FLAGS}")

    # Determine if the processor supports 64bit execution
    EXECUTE_PROCESS(
      COMMAND sysctl hw.cpu64bit_capable
      ERROR_QUIET
      OUTPUT_VARIABLE 64_CMD
      OUTPUT_STRIP_TRAILING_WHITESPACE
      )
    STRING(REGEX REPLACE "^hw.cpu64bit_capable: (.*)" "\\1" 64_BIT "${64_CMD}")
    #  ELSEIF(UNIX)# UNIX must be after APPLE becaus APPLE is UNIX too
  ENDIF()
ENDIF()

IF(NOT CMAKE_COMPILER_IS_CLANG)
  # Compiler options for profiling... but not possible with clang
  OPTION ( PROFILING "Enable/Disable PROFILING" OFF )
  IF(PROFILING)
    ADD_DEFINITIONS(-pg)
    SET(CMAKE_EXE_LINKER_FLAGS "-pg ${CMAKE_EXE_LINKER_FLAGS}")
  ENDIF(PROFILING)
ENDIF(NOT CMAKE_COMPILER_IS_CLANG)

###############################################################################
#####
#####         Choose executable target to compile
#####
###############################################################################
IF(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  MESSAGE(STATUS "Setting build type to 'Release' as none was specified.")
  seT(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  SET_PROPERTY(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
  SET(CMAKE_CONFIGURATION_TYPES ${CMAKE_BUILD_TYPE} )
ENDIF()

############################################################################
#####
#####        MMG (for mesh data structure)
#####
############################################################################
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/modules)

SET(MMG_DIR "" CACHE PATH "Installation directory for mmg")

FIND_PACKAGE(MMG)

IF ( NOT MMG_FOUND )
  MESSAGE ( FATAL_ERROR "ERROR: The installation directory for mmg is required:"
    "(see https://github.com/MmgTools/mmg and  download the branch develop)."
    "If you have already installed Mmg and want to use it, "
    "please set the CMake variable or environment variable MMG_DIR "
    "to your mmg directory.")     # Set flags for building test program
ELSE ( )
  INCLUDE_DIRECTORIES(${MMG_INCLUDE_DIRS})

  MESSAGE(STATUS
    "Compilation with Mmg: ${MMG_LIBRARIES}")
  SET( LIBRARIES ${MMG_LIBRARIES} ${LIBRARIES})

ENDIF ( )

############################################################################
#####
#####         Fortran header: libmmgtypesf.h
#####
############################################################################
ADD_EXECUTABLE(genheader ${PROJECT_SOURCE_DIR}/scripts/genheader.c)
TARGET_LINK_LIBRARIES ( genheader ${LIBRARIES}  )
GENERATE_FORTRAN_HEADER ( hex2tet
  ${H2T_SOURCE_DIR} libhex2tet.h ${MMG_INCLUDE_DIRS}/mmg/common/libmmgtypes.h
  ${H2T_BUILD_DIR} libhex2tetf.h
  )

############################################################################
#####
#####         Scotch
#####
############################################################################
# Find SCOTCH library?
SET(SCOTCH_DIR "" CACHE PATH "Installation directory for scotch")

# add Scotch library?
SET ( USE_SCOTCH "" CACHE STRING "Use SCOTCH tool for renumbering (ON, OFF or <empty>)" )
SET_PROPERTY(CACHE USE_SCOTCH PROPERTY STRINGS "ON" "OFF" "")

IF ( NOT DEFINED USE_SCOTCH OR USE_SCOTCH STREQUAL "" OR USE_SCOTCH MATCHES " +" )
  # Variable is not provided by user
  FIND_PACKAGE(SCOTCH QUIET)

ELSE ()
  IF ( USE_SCOTCH )
    FIND_PACKAGE(SCOTCH)

    IF ( NOT SCOTCH_FOUND )
      MESSAGE ( WARNING "Scotch library not found:"
        "Using scotch reduce the execution time of mmg3d "
        "(see https://gforge.inria.fr/frs/?group_id=248 to download it)."
        "If you have already installed Scotch and want to use it, "
        "please set the CMake variable or environment variable SCOTCH_DIR "
        "to your scotch directory.")
    ENDIF ( )
  ENDIF()
ENDIF()

If ( SCOTCH_FOUND )
  add_definitions(-DUSE_SCOTCH)

  MESSAGE(STATUS
    "Compilation with scotch: ${SCOTCH_LIBRARIES}")
  SET( LIBRARIES  ${LIBRARIES} ${SCOTCH_LIBRARIES})

  INCLUDE_DIRECTORIES(${SCOTCH_INCLUDE_DIRS})

ENDIF()

############################################################################
#####
#####         Elastic
#####
############################################################################
# add SUscElas library?
SET(ELAS_DIR "" CACHE PATH "Installation directory for Elas")
INCLUDE(cmake/modules/FindElas.cmake)

CMAKE_DEPENDENT_OPTION (
  USE_ELAS "Use the Elas library for lagrangian motion option" ON
  "ELAS_FOUND" OFF)

IF( USE_ELAS )
  # Set flags for building test program
  INCLUDE_DIRECTORIES(${ELAS_INCLUDE_DIR})

  SET(CMAKE_REQUIRED_INCLUDES ${ELAS_INCLUDE_DIR})
  SET(CMAKE_REQUIRED_LIBRARIES ${ELAS_LIBRARY})

  SET(CMAKE_C_FLAGS "-DUSE_ELAS ${CMAKE_C_FLAGS}")
  MESSAGE(STATUS
    "Compilation with the Elas library: ${ELAS_LIBRARY} ")
  SET( LIBRARIES ${ELAS_LINK_FLAGS} ${LIBRARIES})
  SET( LIBRARIES ${ELAS_LIBRARY} ${LIBRARIES})
ENDIF()

IF (ELAS_NOTFOUND)
  MESSAGE ( WARNING "Elas is a library to solve the linear elasticity "
    "problem (see https://github.com/SUscTools/Elas to download it). "
    "This library is needed to use the lagrangian motion option. "
    "If you have already installed Elas and want to use it, "
    "please set the CMake variable or environment variable ELAS_DIR "
    "to your Elas directory.")
ENDIF ( )

###############################################################################
#####
#####         Sources and libraries
#####
###############################################################################
# Header files
INCLUDE_DIRECTORIES(${H2T_SOURCE_DIR})
# Source files
FILE(
  GLOB
  library_files
  ${H2T_SOURCE_DIR}/*.c   ${H2T_SOURCE_DIR}/*.h
  )
LIST(REMOVE_ITEM library_files
  ${H2T_SOURCE_DIR}/hex2tet.c)

FILE(
  GLOB
  main_file
  ${H2T_SOURCE_DIR}/hex2tet.c
  )
IF (NOT WIN32)
  FIND_LIBRARY(M_LIB m)
  SET( LIBRARIES ${LIBRARIES} ${M_LIB})
ENDIF()

############################################################################
#####
##### RPATH for MACOSX
#####
############################################################################
SET(CMAKE_MACOSX_RPATH 1)

############################################################################
#####
#####         Compile hex2tet libraries
#####
############################################################################
OPTION ( LIBHEX2TET_STATIC "Compile the hex2tet static library" OFF)
OPTION ( LIBHEX2TET_SHARED "Compile the hex2tet shared library" OFF)

# Compile static library
IF ( LIBHEX2TET_STATIC )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_a STATIC
    "${library_files}" ${PROJECT_NAME} )
ENDIF()

# Compile shared library
IF ( LIBHEX2TET_SHARED )
  ADD_AND_INSTALL_LIBRARY ( lib${PROJECT_NAME}_so SHARED
    "${library_files}" ${PROJECT_NAME} )
ENDIF()

# hex2tet header files needed for library
SET( h2t_headers
  ${H2T_SOURCE_DIR}/libhex2tet.h
  ${H2T_BUILD_DIR}/libhex2tetf.h
  )

# Install header files in /usr/local or equivalent
INSTALL(FILES ${h2t_headers} DESTINATION include/hex2tet COMPONENT headers)

# Copy header files in project directory at every build step
COPY_HEADERS_AND_CREATE_TARGET ( ${H2T_SOURCE_DIR} ${H2T_BUILD_DIR} ${H2T_INCLUDE} hex2tet )


############################################################################
#####
#####         Compile program to test library
#####
############################################################################

IF ( TEST_LIBMMG3D )
  INCLUDE(cmake/testing/h2t_tests.cmake)
ENDIF()

###############################################################################
#####
#####         Compile H2T executable
#####
###############################################################################
ADD_AND_INSTALL_EXECUTABLE ( ${PROJECT_NAME} "${library_files}" ${main_file} )

###############################################################################
#####
#####         Create API Documentation
#####
###############################################################################
FIND_PACKAGE(Doxygen)
IF(DOXYGEN_FOUND)

  # hex2tet Documentation
  CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile.in
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile @ONLY)

  ADD_CUSTOM_TARGET(h2t_doc
    COMMAND ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/Doxyfile
    # COMMAND ${CMAKE_COMMAND} -E chdir
    # ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/latex make
    # COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/latex/refman.pdf
    # ${CMAKE_CURRENT_SOURCE_DIR}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/doc/doxygen/
    COMMENT "Generating hex2tet API documentation with Doxygen.
    Open the doc/doxygen/html/index.html file to see it." VERBATIM
    )

  ADD_CUSTOM_TARGET(doc
    DEPENDS h2t_doc
    COMMENT "Generating hex2tet API documentation with Doxygen.
     Open the doc/doxygen/html/index.html file to see it" VERBATIM
    )
ENDIF ( DOXYGEN_FOUND )

OPTION ( BUILD_TESTING "Enable / Disable tests" OFF )
INCLUDE( ${CMAKE_SOURCE_DIR}/cmake/testing/h2t_tests.cmake )
