#------------------------------------------------------------------------------
# CMakeLists.txt
# Module : G4visualization
# Package: Geant4.src.G4visualization
#
# Intermediate level CMakeLists.txt - just process subdirectories
#
# Generated on : 24/9/2010
#
# $Id: CMakeLists.txt,v 1.1 2010-09-29 19:11:25 bmorgan Exp $
#
#------------------------------------------------------------------------------

add_subdirectory(FukuiRenderer)
add_subdirectory(HepRep)
add_subdirectory(RayTracer)
add_subdirectory(Tree)
add_subdirectory(VRML)
add_subdirectory(XXX)
add_subdirectory(externals)
add_subdirectory(gMocren)
add_subdirectory(management)
add_subdirectory(modeling)

#
# OpenGL is optional depending on user selection
#
if(GEANT4_USE_OPENGL)
    add_subdirectory(OpenGL)
endif()

#
# OpenInventor is optional depending on user selection
#
if(GEANT4_USE_INVENTOR)
    add_subdirectory(OpenInventor)
endif()

