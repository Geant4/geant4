#.rst
# G4ClangFormat
# =============
#
# Functions and helper targets for formatting sources of a target using 
# clang-format

# - Include guard
if(NOT __G4CLANGFORMAT_INCLUDED)
  set(__G4CLANGFORMAT_INCLUDED 1)
else()
  return()
endif()

# Default logging directory
set(G4FORMAT_LOGDIR ${PROJECT_BINARY_DIR}/format CACHE PATH
    "Output folder of files that are formatted")

# Find clang-format (optional)
find_program(CLANG_FORMATTER NAMES
    clang-format-6          # prefer clang-format version 6.0
    clang-format-6.0
    clang-format-mp-6.0     # Apple macports version
    clang-format)

function(exclude_from_format)
  cmake_parse_arguments(ARG "" "" "HEADERS;SOURCES;FILES" ${ARGN})

  function(_add_labels _path)
    foreach(_FILE ${ARGN})
      if(NOT IS_ABSOLUTE ${_FILE})
        set(_FILE ${CMAKE_CURRENT_LIST_DIR}${_path}/${_FILE})
      endif()
      set(_LABELS "EXCLUDE_FORMAT")
      get_source_file_property(_EXISTING_LABELS ${_FILE} LABELS)
      if(_EXISTING_LABELS)
        list(APPEND _LABELS ${_EXISTING_LABELS})
      endif()
      set_source_files_properties(${_FILE} PROPERTIES LABELS "${_LABELS}")
      get_source_file_property(_EXISTING_LABELS ${_FILE} LABELS)
    endforeach()
  endfunction()

  _add_labels("/include" ${ARG_HEADERS})
  _add_labels("/src" ${ARG_SOURCES})
  _add_labels("" ${ARG_FILES})
endfunction()

#------------------------------------------------------------------------------#

function(create_format_directory)
  if(NOT EXISTS ${G4FORMAT_LOGDIR})
      execute_process(COMMAND ${CMAKE_COMMAND}
          -E make_directory ${G4FORMAT_LOGDIR})
  endif()
endfunction()

#------------------------------------------------------------------------------#

function(geant4_format_target)
  if(NOT CLANG_FORMATTER)
    return()
  endif()

  cmake_parse_arguments(G4FMTTARGET "" "NAME" "SOURCES" ${ARGN})
  string(REGEX REPLACE "-format$" "" G4TARGET_NAME "${G4FMTTARGET_NAME}")

  # ensure not empty list
  if("${G4FMTTARGET_SOURCES}" STREQUAL "")
    return()
  endif()

  # create output directory
  create_format_directory()
  set(SOURCES ) # list of sources
  set(FILELOG ) # string of sources (for log)
  foreach(_SOURCE ${G4FMTTARGET_SOURCES})
    # files in binary directory are almost always generated files
    string(FIND "${_SOURCE}" "${PROJECT_BINARY_DIR}" IN_BINARY_DIR)
    # ignore generated files
    if(IN_BINARY_DIR GREATER -1)
      continue()
    endif()
    # check if exists (valid full path)
    if(EXISTS "${_SOURCE}")
      # do nothing -- absolute path was specified
    elseif(EXISTS "${CMAKE_CURRENT_LIST_DIR}/${_SOURCE}")
      # a relative path was specified
      set(_SOURCE ${CMAKE_CURRENT_LIST_DIR}/${_SOURCE})
    else()
      # neither relative nor absolute
      message(STATUS "[format] file not found: '${_SOURCE}'")
      continue()
    endif()
    #string(REPLACE "${PROJECT_SOURCE_DIR}/" "" _SOURCE ${_SOURCE})
    get_source_file_property(SOURCE_LABELS ${_SOURCE} LABELS)
    list(FIND SOURCE_LABELS "EXCLUDE_FORMAT" EXCLUDE_FORMAT)
    if(EXCLUDE_FORMAT GREATER -1)
      continue()
    endif()
    # append to log
    set(FILELOG "${FILELOG}${_SOURCE}\n")
    # add to list to be processed
    list(APPEND SOURCES ${_SOURCE})
  endforeach()
    
  # write format file
  set(FMTLOG_OUT ${G4FORMAT_LOGDIR}/${G4TARGET_NAME}.txt)
  file(WRITE ${FMTLOG_OUT} "${FILELOG}")
  file(RELATIVE_PATH FMTLOG_RELATIVE_OUT ${PROJECT_BINARY_DIR} ${FMTLOG_OUT})
  # get the number of source files
  list(LENGTH SOURCES NUM_SOURCES)
  # if there are too many source files (arguments), command will fail
  # so process no more than 500 at a time
  if(NUM_SOURCES GREATER 500)
    # function to generate indices
    function(generate_indices VAR NUM)
      set(INDICES )
      set(N 0)
      while(N LESS NUM)
        list(APPEND INDICES ${N})
        math(EXPR N "${N}+1")
      endwhile()
      set(${VAR} ${INDICES} PARENT_SCOPE)
    endfunction()
    # generate 500 indices
    generate_indices(INDICES 500)
    set(FORMAT_COMMANDS )
    # while there are still sources in list
    while(NUM_SOURCES GREATER 0)
      # if below 500, re-generate indices
      if(NUM_SOURCES LESS 500)
        generate_indices(INDICES ${NUM_SOURCES})
      endif()
      # get the subset of first 500
      list(GET SOURCES ${INDICES} SOURCES_SUBSET)
      # add the target for block of 500
      list(APPEND FORMAT_COMMANDS COMMAND ${CLANG_FORMATTER} -i ${SOURCES_SUBSET})
      # remove 500 we just processed
      list(REMOVE_AT SOURCES ${INDICES})
      # update the number of sources remaining
      list(LENGTH SOURCES NUM_SOURCES)
    endwhile()
    # add the custom command
    list(LENGTH FORMAT_COMMANDS NUM_COMMANDS)
    if(NUM_COMMANDS GREATER 0)
      add_custom_target(${G4FMTTARGET_NAME}
        ${FORMAT_COMMANDS}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT "Running '${CLANG_FORMATTER}' on source files in '${G4TARGET_NAME}' target (see ${FMTLOG_RELATIVE_OUT})...")
    endif()
  else()
    list(LENGTH SOURCES NUM_SOURCES)
    if(NUM_SOURCES GREATER 0)
      # if less than 500 source files, add all in one target
      add_custom_target(${G4FMTTARGET_NAME}
        COMMAND ${CLANG_FORMATTER} -i ${SOURCES}
        WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
        COMMENT
        "Running '${CLANG_FORMATTER}' on source files in '${G4TARGET_NAME}' target (see ${FMTLOG_RELATIVE_OUT})...")
    endif()
  endif(NUM_SOURCES GREATER 500)
endfunction()
