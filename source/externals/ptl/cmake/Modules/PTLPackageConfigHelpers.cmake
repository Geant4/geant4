# -------------------------------------------------------------------------------------- #
# PTL CMake Package Config installation
#
include(CMakePackageConfigHelpers)

# Generic
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

# Install Tree
set(INCLUDE_INSTALL_DIR "${PTL_INSTALL_INCLUDEDIR}")
set(LIB_INSTALL_DIR "${PTL_INSTALL_LIBDIR}")
set(CMAKE_MODULE_INSTALL_DIR "${PTL_INSTALL_CMAKEDIR}/Modules")

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/${PROJECT_NAME}Config.cmake.in
    ${PROJECT_BINARY_DIR}/installation/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${PTL_INSTALL_CMAKEDIR}
    INSTALL_PREFIX ${PTL_INSTALL_PREFIX}
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR CMAKE_MODULE_INSTALL_DIR)

install(
    FILES ${PROJECT_BINARY_DIR}/installation/${PROJECT_NAME}Config.cmake
          ${PROJECT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
    DESTINATION ${PTL_INSTALL_CMAKEDIR}
    COMPONENT Development)

install(
    FILES ${PROJECT_SOURCE_DIR}/cmake/Modules/FindTBB.cmake
    DESTINATION ${PTL_INSTALL_CMAKEDIR}/Modules
    COMPONENT Development)

# Build Tree
set(INCLUDE_INSTALL_DIR "${PROJECT_SOURCE_DIR}/include")
set(LIB_INSTALL_DIR "${PROJECT_BINARY_DIR}")
set(CMAKE_MODULE_INSTALL_DIR "${PROJECT_SOURCE_DIR}/cmake/Modules")

configure_package_config_file(
    ${PROJECT_SOURCE_DIR}/cmake/Templates/${PROJECT_NAME}Config.cmake.in
    ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
    INSTALL_DESTINATION ${PROJECT_BINARY_DIR}
    INSTALL_PREFIX ${PROJECT_BINARY_DIR}
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR CMAKE_MODULE_INSTALL_DIR)

configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/${PROJECT_NAME}BuildTargets.cmake
               ${PROJECT_BINARY_DIR}/${PROJECT_NAME}Targets.cmake COPYONLY)

# -------------------------------------------------------------------------------------- #
# PTL pkg-config file install NB: The `-std=c++<EPOCH>` flag is currently exported to
# Cflags. Not at all clear how this is supposed to be handled in pkg-config yet. This is
# therefore a "simplest and dumbest" solution
set(PTL_PC_INCLUDEDIR "\${prefix}/${PTL_INSTALL_INCLUDEDIR}")
if(IS_ABSOLUTE "${PTL_INSTALL_INCLUDEDIR}")
    set(PTL_PC_INCLUDEDIR "${PTL_INSTALL_INCLUDEDIR}")
endif()

set(PTL_PC_LIBDIR "\${prefix}/${PTL_INSTALL_LIBDIR}")
if(IS_ABSOLUTE "${PTL_INSTALL_LIBDIR}")
    set(PTL_PC_INCLUDEDIR "${PTL_INSTALL_LIBDIR}")
endif()

if(PTL_USE_TBB)
    # Left as exact until we know version/soversions for TBB better
    set(PTL_PC_TBB_REQUIREMENT "tbb = ${TBB_VERSION}")
endif()

configure_file(${PROJECT_SOURCE_DIR}/cmake/Templates/ptl.pc.in
               ${PROJECT_BINARY_DIR}/installation/ptl.pc @ONLY)

install(
    FILES "${PROJECT_BINARY_DIR}/installation/ptl.pc"
    RENAME "G4ptl.pc"
    DESTINATION ${PTL_INSTALL_LIBDIR}/pkgconfig/
    COMPONENT Development)
