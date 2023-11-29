# - Basic setup for packaging Geant4 using CPack
# Still rather rough and only generally intended to create source packages
# at the moment.
#

#-----------------------------------------------------------------------
# Package up needed system libraries - only for WIN32?
#
# for Intel icc suppress warnings about icc libraries not existing for 19. 
if(CMAKE_CXX_COMPILER_ID MATCHES "Intel")
   set(CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS_NO_WARNINGS True)
endif()

include(G4GitUtilities)
include(InstallRequiredSystemLibraries)

#-----------------------------------------------------------------------
# Copy/Generate common resource files into formats CPack generators like
file(WRITE ${PROJECT_BINARY_DIR}/README.txt "
Geant4
======
A toolkit for the simulation of the passage of particles through matter.
Its areas of application include high energy, nuclear and accelerator physics, as well
as studies in medical and space science.

For more information about how to use Geant4 and its application domains, please visit
the Geant4 website (https://cern.ch/geant4)

There are three main reference papers for Geant4:

- Nuclear Instruments and Methods in Physics Research A 835 (2016) 186-225 
  - http://www.sciencedirect.com/science/article/pii/S0168900216306957
- IEEE Transactions on Nuclear Science 53 No. 1 (2006) 270-278
  - https://ieeexplore.ieee.org/xpls/abs_all.jsp?isnumber=33833&amp;arnumber=1610988&amp;count=33&amp;index=7
- Nuclear Instruments and Methods in Physics Research A 506 (2003) 250-303
  - http://www.sciencedirect.com/science/article/pii/S0168900203013688

This installer has been created using CPack (http://www.cmake.org).
")

configure_file(LICENSE LICENSE.txt COPYONLY) # Suffix must be .txt

#-----------------------------------------------------------------------
# General packaging setup - variables relevant to all package formats
#
set(CPACK_PACKAGE_NAME "geant4")
set(CPACK_PACKAGE_VENDOR "Geant4 Collaboration")
set(CPACK_PACKAGE_CONTACT "Geant4 Collaboration <https://geant4.cern.ch/collaboration/contacts>")
set(CPACK_PACKAGE_VERSION ${Geant4_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${Geant4_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${Geant4_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${Geant4_VERSION_PATCH})
set(CPACK_PACKAGE_DESCRIPTION_FILE "${PROJECT_BINARY_DIR}/README.txt")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Geant4-${Geant4_VERSION_MAJOR}.${Geant4_VERSION_MINOR}")
set(CPACK_PACKAGE_ICON "${PROJECT_SOURCE_DIR}/cmake/Templates/g4_small.bmp")
set(CPACK_PACKAGE_CHECKSUM "SHA256")

set(CPACK_DEBIAN_FILE_NAME "DEB-DEFAULT")
set(CPACK_DEBIAN_PACKAGE_SECTION "science")
set(CPACK_DEBIAN_PACKAGE_SHLIBDEPS TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS TRUE)
set(CPACK_DEBIAN_PACKAGE_GENERATE_SHLIBS_POLICY "=")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "cmake, make, g++")
set(CPACK_DEBIAN_PACKAGE_RECOMMENDS "libexpat1-dev, libxerces-c-dev, zlib1g-dev")
set(CPACK_DEBIAN_PACKAGE_SUGGESTS "")

# - Common Resource files
set(CPACK_RESOURCE_FILE_README  "${PROJECT_BINARY_DIR}/README.txt")
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_BINARY_DIR}/LICENSE.txt")

#-----------------------------------------------------------------------
# Source package settings
#
# - Use same source package naming as GitLab generates for tags
set(CPACK_SOURCE_PACKAGE_FILE_NAME "geant4-v${CPACK_PACKAGE_VERSION}")
set(CPACK_SOURCE_GENERATOR "TGZ;TXZ;ZIP")
set(CPACK_SOURCE_STRIP_FILES "")

# Exclude VCS, temp files, and development-only components from the source package.
set(CPACK_SOURCE_IGNORE_FILES
  "${PROJECT_BINARY_DIR}/"
  ".clang-format"
  ".clang-tidy"
  "README.*"
  "CONTRIBUTING.*"
  "CODING_GUIDELINES.*"
  "RELEASE_MANAGEMENT.rst"
  "SHIFTERS.rst" 
  "/GitUtilities/"
  "/tests/"
  "/ReleaseNotes/development"
  "/benchmarks/"
  "/verification/"
  "/test.*/"
  "/doc.*/"
  "/benchmark.*/"
  "~$"
  "/CVS/"
  "/.svn/"
  "/\\\\\\\\.svn/"
  "/.git/"
  ".gitignore"
  "/.github/"
  "/.gitlab"
  "/\\\\\\\\.git/"
  "\\\\\\\\.swp$"
  "\\\\\\\\.swp$"
  "\\\\.swp"
  "\\\\\\\\.#"
  "/#")

# Add any _configure time_ untracked files. Actual package build will warn/error out
# on any unknown ones. Due to an unknown quirk in CPack, changes to CPACK_SOURCE_IGNORE_FILES
# after inclusion of CPack (e.g. in `CPACK_PROJECT_CONFIG_FILE`) are not taken into account
# when CPack runs.
geant4_git_find_dirty(${PROJECT_SOURCE_DIR} __modded __untracked)
list(APPEND CPACK_SOURCE_IGNORE_FILES ${__untracked})

# Create and inject sources that require generation or otherwise including
configure_file(${PROJECT_SOURCE_DIR}/README.rst _source_extras/README.rst COPYONLY)

if(EXISTS "${PROJECT_SOURCE_DIR}/CONTRIBUTING_GitHub.rst")
  configure_file(${PROJECT_SOURCE_DIR}/CONTRIBUTING_GitHub.rst _source_extras/CONTRIBUTING.rst COPYONLY)
else()
  configure_file(${PROJECT_SOURCE_DIR}/CONTRIBUTING.rst _source_extras/CONTRIBUTING.rst COPYONLY)
endif()



configure_file(
  cmake/Templates/source_package_extras.cmake.in
  source_package_extras.cmake
  @ONLY)
list(APPEND CPACK_INSTALL_SCRIPTS "${CMAKE_CURRENT_BINARY_DIR}/source_package_extras.cmake")

#-----------------------------------------------------------------------
# Binary package common settings
#
set(CPACK_PACKAGE_RELOCATABLE ${CMAKE_INSTALL_IS_NONRELOCATABLE})

if(WIN32)
  set(CPACK_GENERATOR "NSIS;ZIP")
elseif(APPLE)
  set(CPACK_GENERATOR "productbuild;TGZ")
else()
  set(CPACK_GENERATOR "STGZ;TGZ")
endif()

#-----------------------------------------------------------------------
# Finally, generate the CPack per-generator options file and include the
# base CPack configuration.
#
configure_file(
  cmake/Templates/CMakeCPackOptions.cmake.in
  CMakeCPackOptions.cmake
  @ONLY
  )
set(CPACK_PROJECT_CONFIG_FILE ${PROJECT_BINARY_DIR}/CMakeCPackOptions.cmake)
include(CPack)

#-----------------------------------------------------------------------
# Define components and installation types (after CPack included!)
#
cpack_add_install_type(full      DISPLAY_NAME "Full Installation")
cpack_add_install_type(runtime   DISPLAY_NAME "Runtime Installation")
cpack_add_install_type(developer DISPLAY_NAME "Developer Installation")

# - Components for Development
cpack_add_component(Development
  DISPLAY_NAME "Development Components"
  DESCRIPTION "Install all files needed for developing Geant4 applications (headers, makefiles, etc.)"
  INSTALL_TYPES developer full
  )

# - Components for Runtime
cpack_add_component(Runtime
  DISPLAY_NAME "Geant4 runtime Libraries"
  DESCRIPTION "Install all Geant4 libraries"
  INSTALL_TYPES runtime developer full
  )

# - Components for Examples
cpack_add_component(Examples
  DISPLAY_NAME "Application Examples"
  DESCRIPTION "Install all Geant4 examples"
  INSTALL_TYPES full developer
  )

# - Components for Data
cpack_add_component(Data
  DISPLAY_NAME "Geant4 Data Files"
  DESCRIPTION "Install all Geant4 data files"
  INSTALL_TYPES full
  )

