###############################################################################
#####
#####         Generation of a fortran header file
#####
###############################################################################

MACRO ( GENERATE_FORTRAN_HEADER name
    in_dir in_file out_dir out_file
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT ${out_dir}/${out_file}
    COMMAND genheader ${out_dir}/${out_file} ${in_dir}/${in_file}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS genheader ${in_dir}/${in_file}
    ${PROJECT_SOURCE_DIR}/scripts/genfort.pl
    COMMENT "Generating Fortran header for ${name}"
    )

  ADD_CUSTOM_TARGET (
    ${name}_fortran_header
    ALL
    DEPENDS ${out_dir}/${out_file}
    )

ENDMACRO ( )

###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####
###############################################################################
MACRO ( COPY_HEADER
    in_dir in_file out_dir out_file
    file_dependencies
    target_name
    )
  # Wrap add_custom_command into add_custom target to remove dpendencies from
  # the custom command and thus allow parallel build.
  ADD_CUSTOM_COMMAND (
    OUTPUT  ${out_dir}/${out_file}
    COMMAND ${CMAKE_COMMAND} -E copy  ${in_dir}/${in_file} ${out_dir}/${out_file}
    WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    DEPENDS ${file_dependencies} ${in_dir}/${in_file}
    COMMENT "Copying ${in_dir}/${in_file} in ${out_dir}/${out_file}"
    )

  ADD_CUSTOM_TARGET ( ${target_name} ALL
    DEPENDS ${out_dir}/${out_file} )
  ADD_DEPENDENCIES ( ${target_name} ${file_dependencies} )

ENDMACRO ( )


###############################################################################
#####
#####         Copy an automatically generated header file to another place
#####         and create the associated target
#####
###############################################################################
MACRO ( COPY_HEADERS_AND_CREATE_TARGET
    source_dir binary_dir include_dir target_identifier )

  ADD_CUSTOM_TARGET(${target_identifier}_header ALL
    DEPENDS
    ${source_dir}/lib${target_identifier}.h )

  COPY_HEADER (
    ${source_dir} lib${target_identifier}.h
    ${include_dir} lib${target_identifier}.h
    ${target_identifier}_header copy_lib${target_identifier}
    )

  COPY_HEADER (
    ${binary_dir} lib${target_identifier}f.h
    ${include_dir} lib${target_identifier}f.h
    ${target_identifier}_fortran_header copy_lib${target_identifier}f
    )

  ADD_CUSTOM_TARGET(copy_${target_identifier}_headers ALL
    DEPENDS
    copy_lib${target_identifier}f
    copy_lib${target_identifier} )


ENDMACRO ( )

###############################################################################
#####
#####         Add a library to build and needed include dir, set its
#####         properties, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_LIBRARY
    target_name target_type sources output_name )

  ADD_LIBRARY ( ${target_name} ${target_type} ${sources} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PUBLIC
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} PUBLIC
      ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
  ENDIF ( )

  SET_TARGET_PROPERTIES ( ${target_name}
    PROPERTIES OUTPUT_NAME ${output_name} )

  SET_PROPERTY(TARGET ${target_name} PROPERTY C_STANDARD 99)

  TARGET_LINK_LIBRARIES ( ${target_name} ${LIBRARIES} )

  INSTALL ( TARGETS ${target_name}
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    COMPONENT lib)

ENDMACRO ( )

###############################################################################
#####
#####         Add an executable to build and needed include dir, set its
#####         postfix, add link dependencies and the install rule
#####
###############################################################################

MACRO ( ADD_AND_INSTALL_EXECUTABLE
    exec_name lib_files main_file )

  IF ( NOT TARGET lib${exec_name}_a AND NOT TARGET lib${exec_name}_so )
    ADD_EXECUTABLE ( ${exec_name} ${lib_files} ${main_file} )
  ELSE ( )
    ADD_EXECUTABLE ( ${exec_name} ${main_file})

    SET_PROPERTY(TARGET ${exec_name} PROPERTY C_STANDARD 99)

    IF ( NOT TARGET lib${exec_name}_a )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_so)
    ELSE ( )
      TARGET_LINK_LIBRARIES(${exec_name} lib${exec_name}_a)
    ENDIF ( )

  ENDIF ( )

  IF ( WIN32 AND NOT MINGW AND USE_SCOTCH )
    my_add_link_flags ( ${exec_name} "/SAFESEH:NO")
  ENDIF ( )

 IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
   INCLUDE_DIRECTORIES ( ${exec_name} PUBLIC
     ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
 ELSE ( )
   TARGET_INCLUDE_DIRECTORIES ( ${exec_name} PUBLIC
     ${COMMON_BINARY_DIR} ${COMMON_SOURCE_DIR} ${PROJECT_BINARY_DIR}/include )
 ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${exec_name} ${LIBRARIES}  )

  INSTALL(TARGETS ${exec_name} RUNTIME DESTINATION bin COMPONENT appli)

  ADD_TARGET_POSTFIX(${exec_name})

ENDMACRO ( )


###############################################################################
#####
#####         Add a target postfix depending on the build type
#####
###############################################################################

MACRO ( ADD_TARGET_POSTFIX target_name )
  # in debug mode we name the executable mmgs_debug
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES DEBUG_POSTFIX _debug)
  # in Release mode we name the executable mmgs_O3
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELEASE_POSTFIX _O3)
  # in RelWithDebInfo mode we name the executable mmgs_O3d
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES RELWITHDEBINFO_POSTFIX _O3d)
  # in MinSizeRel mode we name the executable mmgs_O3
  SET_TARGET_PROPERTIES(${target_name} PROPERTIES MINSIZEREL_POSTFIX _Os)
ENDMACRO ( )

###############################################################################
#####
#####         Add Executable that must be tested by ci
#####
###############################################################################

MACRO ( ADD_EXEC_TO_CI_TESTS exec_name list_name )

  IF(${CMAKE_BUILD_TYPE} MATCHES "Debug")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_debug)
    SET(BUILDNAME ${BUILDNAME}_debug CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "Release")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_O3)
    SET(BUILDNAME ${BUILDNAME}_O3 CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "RelWithDebInfo")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_O3d)
    SET(BUILDNAME ${BUILDNAME}_O3d CACHE STRING "build name variable")
  ELSEIF(${CMAKE_BUILD_TYPE} MATCHES "MinSizeRel")
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name}_Os)
    SET(BUILDNAME ${BUILDNAME}_Os CACHE STRING "build name variable")
  ELSE()
    SET(EXECUT ${EXECUTABLE_OUTPUT_PATH}/${exec_name})
    SET(BUILDNAME ${BUILDNAME} CACHE STRING "build name variable")
  ENDIF()

  SET ( ${list_name} ${EXECUT} )

ENDMACRO ( )


###############################################################################
#####
#####         Add a library test
#####
###############################################################################

MACRO ( ADD_LIBRARY_TEST target_name main_path target_dependency lib_name )
  ADD_EXECUTABLE ( ${target_name} ${main_path} )
  ADD_DEPENDENCIES( ${target_name} ${target_dependency} )

  IF ( CMAKE_VERSION VERSION_LESS 2.8.12 )
    INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${PROJECT_BINARY_DIR}/include )
  ELSE ( )
    TARGET_INCLUDE_DIRECTORIES ( ${target_name} PUBLIC ${PROJECT_BINARY_DIR}/include )
  ENDIF ( )

  IF ( WIN32 AND ((NOT MINGW) AND USE_SCOTCH) )
    MY_ADD_LINK_FLAGS ( ${target_name} "/SAFESEH:NO" )
  ENDIF ( )

  TARGET_LINK_LIBRARIES ( ${target_name}  ${lib_name} )
  INSTALL(TARGETS ${target_name} RUNTIME DESTINATION bin COMPONENT appli )

ENDMACRO ( )
