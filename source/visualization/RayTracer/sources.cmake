#------------------------------------------------------------------------------
# sources.cmake
# Module : G4RayTracer
# Package: Geant4.src.G4visualization.G4RayTracer
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 88190 2015-02-02 17:24:54Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})
include_directories(${USOLIDS_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/detector/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/digits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/digits_hits/hits/include)
include_directories(${CMAKE_SOURCE_DIR}/source/event/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/navigation/include)
include_directories(${CMAKE_SOURCE_DIR}/source/geometry/volumes/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/graphics_reps/include)
include_directories(${CMAKE_SOURCE_DIR}/source/intercoms/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/bosons/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/cuts/include)
include_directories(${CMAKE_SOURCE_DIR}/source/processes/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/run/include)
include_directories(${CMAKE_SOURCE_DIR}/source/track/include)
include_directories(${CMAKE_SOURCE_DIR}/source/tracking/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/visualization/modeling/include)


#
# Module has optional sources
#
# List those always built
set(G4VIS_RAYTRACER_MODULE_HEADERS
    G4RTJpeg.hh
    G4RTJpegCoder.hh
    G4RTJpegCoderTables.hh
    G4RTJpegMaker.hh
    G4RTMessenger.hh
    G4RTOutBitStream.hh
    G4RTPrimaryGeneratorAction.hh
    G4RTRun.hh
    G4RTRunAction.hh
    G4RTSimpleScanner.hh
    G4RTSteppingAction.hh
    G4RTTrackingAction.hh
    G4RTWorkerInitialization.hh
    G4RayTracer.hh
    G4RayTracerFeatures.hh
    G4RayTracerSceneHandler.hh
    G4RayTracerViewer.hh
    G4RayTrajectory.hh
    G4RayTrajectoryPoint.hh
    G4TheMTRayTracer.hh
    G4TheRayTracer.hh
    G4VFigureFileMaker.hh
    G4VRTScanner.hh)

set(G4VIS_RAYTRACER_MODULE_SOURCES
    G4RTJpegCoder.cc
    G4RTJpegMaker.cc
    G4RTMessenger.cc
    G4RTOutBitStream.cc
    G4RTPrimaryGeneratorAction.cc
    G4RTRun.cc
    G4RTRunAction.cc
    G4RTSimpleScanner.cc
    G4RTSteppingAction.cc
    G4RTTrackingAction.cc
    G4RTWorkerInitialization.cc
    G4RayTracer.cc
    G4RayTracerSceneHandler.cc
    G4RayTracerViewer.cc
    G4RayTrajectory.cc
    G4RayTrajectoryPoint.cc
    G4TheMTRayTracer.cc
    G4TheRayTracer.cc
    G4VRTScanner.cc)

set(G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES )

#
# X11 RayTracer only if selected
#
if(GEANT4_USE_RAYTRACER_X11)
    list(APPEND G4VIS_RAYTRACER_MODULE_HEADERS
        G4RayTracerX.hh
        G4RayTracerXViewer.hh
        G4RTXScanner.hh)

    list(APPEND G4VIS_RAYTRACER_MODULE_SOURCES
        G4RayTracerX.cc
        G4RayTracerXViewer.cc
        G4RTXScanner.cc)

    #
    # Add source properties and additional LINK_LIBRARIES here
    #
    # Need X11 includes!
    include_directories(${X11_INCLUDE_DIR})

    # Must use G4VIS_BUILD_RAYTRACERX_DRIVER define
    GEANT4_ADD_COMPILE_DEFINITIONS(SOURCES ${G4VIS_RAYTRACER_MODULE_SOURCES}
        COMPILE_DEFINITIONS G4VIS_BUILD_RAYTRACERX_DRIVER)

    # The X11 Libraries
    list(APPEND G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES ${X11_LIBRARIES})
endif()
    

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4RayTracer
    HEADERS
        ${G4VIS_RAYTRACER_MODULE_HEADERS}
    SOURCES
        ${G4VIS_RAYTRACER_MODULE_SOURCES}
    GRANULAR_DEPENDENCIES
        G4bosons
        G4cuts
        G4detector
        G4digits
        G4event
        G4geometrymng
        G4globman
        G4graphics_reps
        G4hits
        G4intercoms
        G4materials
        G4modeling
        G4navigation
        G4partman
        G4procman
        G4track
        G4tracking
        G4vis_management
        G4volumes
    GLOBAL_DEPENDENCIES
        G4digits_hits
        G4event
        G4geometry
        G4global
        G4graphics_reps
        G4intercoms
        G4materials
        G4modeling
        G4particles
        G4processes
        G4track
        G4tracking
        G4vis_management
    LINK_LIBRARIES
        ${G4VIS_RAYTRACER_MODULE_LINK_LIBRARIES}
)

# List any source specific properties here

