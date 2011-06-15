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
set(CPACK_PACKAGE_VENDOR "Geant4 Collaboration")

set(CPACK_PACKAGE_VERSION ${${PROJECT_NAME}_VERSION})
set(CPACK_PACKAGE_VERSION_MAJOR ${${PROJECT_NAME}_VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${${PROJECT_NAME}_VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${${PROJECT_NAME}_VERSION_PATCH})

# - Resource files
set(CPACK_RESOURCE_FILE_LICENSE "${PROJECT_SOURCE_DIR}/LICENSE")


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

#----------------------------------------------------------------------------
# Binary package setup
# TODO...


#----------------------------------------------------------------------------
# Finally, include the base CPack configuration

# FINAL step - include base CPack
include(CPack)
