# - Basic setup for packaging Geant4 using CPack
# Still rather rough and only generally intended to create source packages
# at the moment.
# 

#----------------------------------------------------------------------------
# Package up needed system libraries - only for WIN32?
#
include(InstallRequiredSystemLibraries)


#----------------------------------------------------------------------------
# General packaging setup - variable relavant to all package formats
#
set(CPACK_PACKAGE_DESCRIPTION "Geant4 Toolkit")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Geant4 Toolkit")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_BINARY_DIR}/README.txt")
set(CPACK_PACKAGE_VENDOR "Geant4 Collaboration")

set(CPACK_PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})

# - Resource files
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_BINARY_DIR}/LICENSE.txt")
set(CPACK_RESOURCE_FILE_README "${CMAKE_BINARY_DIR}/README.txt")


configure_file(LICENSE LICENSE.txt COPYONLY) # Suffix must be .txt 
file(WRITE ${CMAKE_BINARY_DIR}/README.txt "
Geant4

A toolkit for the simulation of the passage of particles through matter. 

Its areas of application include high energy, nuclear and accelerator physics, as well as studies in medical and space science (http://www.geant4.org).
The two main reference papers for Geant4 are published in Nuclear Instruments and Methods in Physics Research A 506 (2003) 250-303, and IEEE Transactions on Nuclear Science 53 No. 1 (2006) 270-278.

This installer has been created using CPack (http://www.cmake.org). 
")

#----------------------------------------------------------------------------
# Source package settings
# Exclude VCS and standard temporary files from the source package.
# Not totally perfected yet!
set(CPACK_SOURCE_IGNORE_FILES 
    ${PROJECT_BINARY_DIR}
    ${PROJECT_SOURCE_DIR}/tests
    "~$"
    "/CVS/"
    "/.svn/"
    "/\\\\\\\\.svn/"
    "/.git/"
    "/\\\\\\\\.git/"
    "\\\\\\\\.swp$"
    "\\\\\\\\.swp$"
    "\\\\.swp"
    "\\\\\\\\.#"
    "/#"
)

# - Ensure all standard source packages are built
set(CPACK_SOURCE_GENERATOR "TGZ;TBZ2;ZIP")
set(CPACK_SOURCE_STRIP_FILES "")

#----------------------------------------------------------------------------
# Binary package setup
#
set(CPACK_PACKAGE_RELOCATABLE True)
set(CPACK_PACKAGE_INSTALL_DIRECTORY "Geant4 ${Geant4_VERSION_MAJOR}.${Geant4_VERSION_MINOR}")
if(CMAKE_BUILD_TYPE STREQUAL Release)
  set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${Geant4_VERSION}-${CMAKE_SYSTEM_NAME}")
else()
  set(CPACK_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${Geant4_VERSION}-${CMAKE_SYSTEM_NAME}-${CMAKE_BUILD_TYPE}")
endif()

if(WIN32)
  set(CPACK_GENERATOR "NSIS;ZIP")  
elseif(APPLE)
  set(CPACK_GENERATOR "PackageMaker;TGZ")
else()
  set(CPACK_GENERATOR "STGZ;TGZ")
endif()

#----------------------------------------------------------------------------
# Finally, generate the CPack per-generator options file and include the
# base CPack configuration.
#
configure_file(cmake/Templates/CMakeCPackOptions.cmake.in CMakeCPackOptions.cmake @ONLY)
set(CPACK_PROJECT_CONFIG_FILE ${CMAKE_BINARY_DIR}/CMakeCPackOptions.cmake)

# FINAL step - include base CPack
include(CPack)

#----------------------------------------------------------------------------
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

