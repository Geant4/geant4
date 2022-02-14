# - G4RayTracer category build

# Add allocation export symbol for the RayTracer module
# NB: we don't add this at the category level because RayTracer is
# the only vis driver that uses exports at the moment. We can move the
# add_definitions higher later if needed.
add_definitions(-DG4VIS_ALLOC_EXPORT)

geant4_global_library_target(COMPONENTS sources.cmake)
