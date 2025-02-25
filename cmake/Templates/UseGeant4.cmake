# - Use file for Geant4
# Add Module directory so that examples can use internal "FindXXX"
# modules. Appended to minimize any conflict with consumer project
# settings
# DEPRECATED: Only needed by certain examples
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/Modules)
