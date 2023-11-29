# - Use file for Geant4
#
# Optional inclusion of internal Use file. This file can contain 
# variables, functions and macros for strict internal use in Geant4, 
# such as building and running validation tests.
#
include(${CMAKE_CURRENT_LIST_DIR}/UseGeant4_internal.cmake OPTIONAL)

# Add Module directory so that examples can use internal "FindXXX"
# modules. Appended to minimize any conflict with consumer project
# settings
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Modules)
