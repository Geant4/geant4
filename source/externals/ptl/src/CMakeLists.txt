# -------------------------------------------------------------------------------------- #
# General
# -------------------------------------------------------------------------------------- #
# Locate sources and headers for this project - headers are included so they will show up
# in IDEs
file(GLOB_RECURSE ptl_headers ${PROJECT_SOURCE_DIR}/include/PTL/*.hh
     ${PROJECT_SOURCE_DIR}/include/PTL/*.icc)
file(GLOB_RECURSE ptl_sources ${CMAKE_CURRENT_LIST_DIR}/*.cc)

# -------------------------------------------------------------------------------------- #
# Config, Version
# -------------------------------------------------------------------------------------- #

configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Config.hh.in
               ${CMAKE_CURRENT_BINARY_DIR}/PTL/Config.hh @ONLY)

list(APPEND ptl_headers ${CMAKE_CURRENT_BINARY_DIR}/PTL/Config.hh)

configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/Version.hh.in
               ${CMAKE_CURRENT_BINARY_DIR}/PTL/Version.hh @ONLY)

list(APPEND ptl_headers ${CMAKE_CURRENT_BINARY_DIR}/PTL/Version.hh)

# -------------------------------------------------------------------------------------- #
# PTL Library
# -------------------------------------------------------------------------------------- #

if(BUILD_OBJECT_LIBS)

    ptl_build_library(
        TYPE OBJECT
        TARGET_NAME ptl-object
        OUTPUT_NAME G4ptl-object
        SOURCES ${ptl_headers} ${ptl_sources})

    target_link_libraries(ptl-object PUBLIC Threads::Threads)
    if(PTL_USE_TBB)
        target_link_libraries(ptl-object PUBLIC TBB::tbb)
    endif()

    target_include_directories(
        ptl-object PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/source>
                          $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>)
endif()

if(BUILD_SHARED_LIBS)

    ptl_build_library(
        TYPE SHARED
        TARGET_NAME ptl-shared
        OUTPUT_NAME G4ptl
        SOURCES ${ptl_headers} ${ptl_sources})

    target_link_libraries(ptl-shared PUBLIC Threads::Threads)
    if(PTL_USE_TBB)
        target_link_libraries(ptl-shared PUBLIC TBB::tbb)
    endif()

    target_compile_definitions(ptl-shared PUBLIC PTL_BUILD_DLL)

    target_include_directories(
        ptl-shared
        PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
               $<INSTALL_INTERFACE:${PTL_INSTALL_INCLUDEDIR}>)

    if(PTL_USE_SANITIZER AND NOT ("${PTL_SANITIZER_TYPE}" STREQUAL "leak"))
        string(REPLACE " " ";" _sanitize_args "${ptl_sanitize_args}")
        set_target_properties(ptl-shared PROPERTIES INTERFACE_LINK_OPTIONS
                                                    "${_sanitize_args}")
    endif()

    export(
        TARGETS ptl-shared
        NAMESPACE PTL::
        FILE ${PROJECT_BINARY_DIR}/ptl-shared.cmake)
endif()

if(BUILD_STATIC_LIBS)

    ptl_build_library(
        TYPE STATIC
        TARGET_NAME ptl-static
        OUTPUT_NAME $<IF:$<PLATFORM_ID:Windows>,G4ptl-static,G4ptl>
        SOURCES ${ptl_headers} ${ptl_sources})

    target_link_libraries(ptl-static PUBLIC Threads::Threads)
    if(PTL_USE_TBB)
        target_link_libraries(ptl-static PUBLIC TBB::tbb)
    endif()

    target_include_directories(
        ptl-static
        PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
               $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
               $<INSTALL_INTERFACE:${PTL_INSTALL_INCLUDEDIR}>)

    export(
        TARGETS ptl-static
        NAMESPACE PTL::
        FILE ${PROJECT_BINARY_DIR}/ptl-static.cmake)
endif()

# -------------------------------------------------------------------------------------- #
# Installation
# -------------------------------------------------------------------------------------- #

if(PTL_INSTALL_CONFIG)
    # install export
    install(
        EXPORT ${PROJECT_NAME}Targets
        NAMESPACE PTL::
        DESTINATION ${PTL_INSTALL_CMAKEDIR}
        COMPONENT Development)
endif()

if(PTL_INSTALL_HEADERS)
    # headers
    install(
        FILES ${ptl_headers}
        DESTINATION ${PTL_INSTALL_INCLUDEDIR}/PTL
        COMPONENT Development)
endif()
