#
# Create a "make doc" target using Doxygen.
#
# Prototype:
#
# GENERATE_DOCUMENTATION(doxygen_config_file)
#
# Parameters:
#
# * doxygen_config_file: Doxygen configuration file (must in the root of the source
#   directory)

include(MacroUtilities)

# -------------------------------------------------------------------------------------- #

# if BUILD_DOXYGEN_DOCS = ON, we want to build docs quietly else, don't build quietly
add_option(${PROJECT_NAME}_DOCS_QUIET "Suppress standard output when making the docs" ON)
mark_as_advanced(${PROJECT_NAME}_DOCS_QUIET)

if(${PROJECT_NAME}_DOCS_QUIET)
    set(DOXYGEN_QUIET YES)
else()
    set(DOXYGEN_QUIET NO)
endif()

# GraphViz dot program is used to build call graphs, caller graphs, class graphs
find_program(GRAPHVIZ_DOT_PATH dot)
mark_as_advanced(GRAPHVIZ_DOT_PATH)

if("${GRAPHVIZ_DOT_PATH}" STREQUAL "GRAPHVIZ_DOT_PATH-NOTFOUND")
    set(DOXYGEN_DOT_FOUND NO)
    set(GRAPHVIZ_DOT_PATH "")
else()
    set(DOXYGEN_DOT_FOUND YES)
    set(GRAPHVIZ_DOT_PATH ${GRAPHVIZ_DOT_PATH})
endif()

# available Doxygen doc formats
set(AVAILABLE_DOXYGEN_DOC_FORMATS HTML LATEX MAN XML RTF)
# we want HTML, LATEX, and MAN to be default
set(_default_on "MAN")
foreach(_doc_format ${AVAILABLE_DOXYGEN_DOC_FORMATS})
    # find if we want it on
    string(REGEX MATCH "${_doc_format}" SET_TO_ON "${_default_on}")
    # if doc format is MAN and it is not a UNIX machine --> turn off
    if("${_doc_format}" STREQUAL "MAN" AND NOT UNIX)
        set(SET_TO_ON "")
    endif()
    # set ON/OFF
    if("${SET_TO_ON}" STREQUAL "")
        set(_default "OFF")
    else()
        set(_default "ON")
    endif()
    # add option
    add_option(${PROJECT_NAME}_${_doc_format}_DOCS
               "Build documentation with ${_doc_format} format" ${_default})
    mark_as_advanced(${PROJECT_NAME}_${_doc_format}_DOCS)
endforeach()

# loop over doc formats and set GENERATE_DOXYGEN_${_doc_format}_DOC to YES/NO
# GENERATE_DOXYGEN_${_doc_format}_DOC is used in configure_file @ONLY
foreach(_doc_format ${AVAILABLE_DOXYGEN_DOC_FORMATS})
    if(${PROJECT_NAME}_${_doc_format}_DOCS)
        set(GENERATE_DOXYGEN_${_doc_format}_DOC YES)
    else()
        set(GENERATE_DOXYGEN_${_doc_format}_DOC NO)
    endif()
endforeach()

if(DOXYGEN_DOT_FOUND)
    set(DOXYGEN_DOT_GRAPH_TYPES CLASS CALL CALLER)
    # options to turn generation of class, call, and caller graphs
    foreach(_graph_type ${DOXYGEN_DOT_GRAPH_TYPES})
        # create CMake doc string
        string(TOLOWER _graph_type_desc ${_graph_type})
        # add option
        option(${PROJECT_NAME}_${_graph_type}_GRAPH "${_message}" ON)
        mark_as_advanced(${PROJECT_NAME}_${_graph_type}_GRAPH)
        # set GENERATE_DOXYGEN_${_graph_type}_GRAPH to YES/NO
        # GENERATE_DOXYGEN_${_graph_type}_GRAPH is used in configure_file @ONLY
        if(${PROJECT_NAME}_${_graph_type}_GRAPH)
            set(GENERATE_DOXYGEN_${_graph_type}_GRAPH YES)
        else()
            set(GENERATE_DOXYGEN_${_graph_type}_GRAPH NO)
        endif()
    endforeach()
endif()

# get the documentation include directories
get_property(BUILDTREE_DIRS GLOBAL PROPERTY PTL_DOCUMENTATION_DIRS)

if("${BUILDTREE_DIRS}" STREQUAL "")
    message(FATAL_ERROR "Property PTL_DOCUMENTATION_DIRS is empty")
endif()

list(REMOVE_DUPLICATES BUILDTREE_DIRS)

set(SOURCE_FILES)
set(EXTENSIONS
    h
    hh
    hpp
    c
    cc
    cpp
    icc
    tcc
    py)
foreach(_DIR ${BUILDTREE_DIRS})
    foreach(_EXT ${EXTENSIONS})
        file(GLOB _FILES "${_DIR}/*.${_EXT}")
        if(_FILES)
            list(APPEND SOURCE_FILES ${_FILES})
        endif()
        unset(_FILES CACHE)
    endforeach()
endforeach()

# Doxyfiles was spaces not semi-colon separated lists
string(REPLACE ";" " " BUILDTREE_DIRS "${BUILDTREE_DIRS}")
string(REPLACE ";" " " SOURCE_FILES "${SOURCE_FILES}")

# -----------------------------------------------------------------------

if(XCODE)
    set(GENERATE_DOCSET_IF_XCODE YES)
else()
    set(GENERATE_DOCSET_IF_XCODE NO)
endif()

# -----------------------------------------------------------------------

configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Doxyfile.in
               ${PROJECT_BINARY_DIR}/doc/Doxyfile.${PROJECT_NAME} @ONLY)

if(${PROJECT_NAME}_HTML_DOCS)
    file(WRITE ${PROJECT_BINARY_DIR}/doc/${PROJECT_NAME}_Documentation.html
         "<meta http-equiv=\"refresh\" content=\"1;url=html/index.html\">")
endif()

# -------------------------------------------------------------------------------------- #
# Macro to generate documentation from:
# http://www.cmake.org/pipermail/cmake/2007-May/014174.html
macro(GENERATE_DOCUMENTATION DOXYGEN_CONFIG_FILE)

    find_package(Doxygen)
    if(NOT Doxygen_FOUND)
        message(
            STATUS
                "Doxygen executable cannot be found. Disable ${PROJECT_NAME}_DOXYGEN_DOCS"
            )
        return()
    endif()
    set(DOXYFILE_FOUND false)

    if(EXISTS ${PROJECT_BINARY_DIR}/doc/${DOXYGEN_CONFIG_FILE})
        set(DOXYFILE_FOUND true)
    else()
        message(
            STATUS
                "Doxygen config file was not found at ${PROJECT_BINARY_DIR}/doc/${DOXYGEN_CONFIG_FILE}"
            )
    endif()

    if(DOXYGEN_FOUND)
        if(DOXYFILE_FOUND)
            # Add target
            add_custom_target(docs ${DOXYGEN_EXECUTABLE}
                                   "${PROJECT_BINARY_DIR}/doc/${DOXYGEN_CONFIG_FILE}")

            install(
                DIRECTORY ${PROJECT_BINARY_DIR}/doc/man/
                DESTINATION ${CMAKE_INSTALL_MANDIR}
                OPTIONAL
                COMPONENT documentation)

            install(
                DIRECTORY ${PROJECT_BINARY_DIR}/doc/html
                DESTINATION ${CMAKE_INSTALL_DOCDIR}
                OPTIONAL
                COMPONENT documentation)

            install(
                DIRECTORY ${PROJECT_BINARY_DIR}/doc/latex
                DESTINATION ${CMAKE_INSTALL_DOCDIR}
                OPTIONAL
                COMPONENT documentation)

            install(
                FILES ${PROJECT_BINARY_DIR}/doc/${PROJECT_NAME}_Documentation.html
                DESTINATION ${CMAKE_INSTALL_DOCDIR}
                OPTIONAL
                COMPONENT documentation)

        else()
            message(
                STATUS
                    "Doxygen configuration file not found - Documentation will not be generated"
                )
        endif()
    else()
        message(STATUS "Doxygen not found - Documentation will not be generated")
    endif()

endmacro()

# -------------------------------------------------------------------------------------- #
# Macro to generate PDF manual from LaTeX using pdflatex assumes manual is in
# ${CMAKE_SOURCE_DIR}/doc
macro(GENERATE_MANUAL MANUAL_TEX MANUAL_BUILD_PATH EXTRA_FILES_TO_COPY)

    find_program(PDFLATEX pdflatex)

    if(PDFLATEX AND NOT "${PDFLATEX}" STREQUAL "PDFLATEX-NOTFOUND")
        # name with no path is given
        set(MANUAL_NAME ${MANUAL_TEX})
        # set to full path
        set(MANUAL_BUILD_PATH ${CMAKE_BINARY_DIR}/${MANUAL_BUILD_PATH})

        if(NOT EXISTS ${CMAKE_SOURCE_DIR}/doc/${MANUAL_TEX})
            message(
                FATAL_ERROR
                    "LaTeX of manual for ${PROJECT_NAME} is not in ${CMAKE_SOURCE_DIR}/doc"
                )
        endif()

        configure_file(${CMAKE_SOURCE_DIR}/doc/${MANUAL_TEX}
                       ${MANUAL_BUILD_PATH}/${MANUAL_NAME} COPYONLY)

        foreach(_file ${EXTRA_FILES_TO_COPY})
            configure_file(${CMAKE_SOURCE_DIR}/doc/${_file} ${MANUAL_BUILD_PATH}/${_file}
                           COPYONLY)
        endforeach()

        add_custom_target(
            man
            ${PDFLATEX} "${MANUAL_NAME}"
            WORKING_DIRECTORY ${MANUAL_BUILD_PATH})
    endif()
endmacro()
