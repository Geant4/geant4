#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4event
# Package: Geant4.src.G4event
#------------------------------------------------------------------------------

# Add allocation export symbol for the event category
add_definitions(-DG4EVENT_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)


