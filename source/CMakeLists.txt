#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : 
# Package: Geant4.src.
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

add_subdirectory(analysis)
add_subdirectory(digits_hits)
add_subdirectory(error_propagation)
add_subdirectory(event)
add_subdirectory(externals)

# Target G4g3tog4 is optional based on user selection...
if(GEANT4_USE_G3TOG4)
    add_subdirectory(g3tog4)
endif()

add_subdirectory(geometry)
add_subdirectory(global)
add_subdirectory(graphics_reps)
add_subdirectory(intercoms)
add_subdirectory(interfaces)
add_subdirectory(materials)
add_subdirectory(parameterisations)
add_subdirectory(particles)
add_subdirectory(persistency)
add_subdirectory(physics_lists)
add_subdirectory(processes)
add_subdirectory(readout)
add_subdirectory(run)
add_subdirectory(track)
add_subdirectory(tracking)
add_subdirectory(visualization)

