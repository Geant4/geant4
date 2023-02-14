cmake_minimum_required(VERSION 3.16...3.21)
project(exampleVecGeomNav)

option(WITH_ROOT "WITH ROOT SUPPORT" OFF)
option(BUILTIN_G4VECGEOMNAV "Build G4VecGeomNav from source (requires VecGeom)" ON)

find_package(Geant4 REQUIRED gdml)
find_package(VecGeom REQUIRED)
if(WITH_ROOT)
  find_package(ROOT REQUIRED)
  find_package(TGeo2VecGeom REQUIRED)
endif()

#-------------------------------------------------------------------
## Option to use built-in G4VecGeomNav
################################################################################
# Minimum version of G4VecGeomNav we need.
set(G4VecGeomNav_VERSION "") # 2021.06.18") # 0.1")

if(BUILTIN_G4VECGEOMNAV)
  include(FetchContent)

  FetchContent_Declare(
    g4vecgeomnav
    GIT_REPOSITORY https://gitlab.cern.ch/VecGeom/g4vecgeomnav.git
    GIT_TAG master
  )

  FetchContent_GetProperties(g4vecgeomnav)

  # G4VecGeomNav builds a shared lib by default. As this links to Geant4 libs, if Geant4 is 
  # static-only, then we should also build G4VecGeomNav statically because:
  # - Geant4 static libs are not guaranteed to be PIC
  # - Geant4 symbols will be duplicated in g4vecgeomnav and the exampleVecGeomNac program
  if(Geant4_static_FOUND AND (NOT Geant4_shared_FOUND))
    set(BUILD_STATIC_LIBS ON)
  endif()

  if(NOT g4vecgeomnav_POPULATED)
    FetchContent_Populate(g4vecgeomnav)
    message(STATUS "${g4vecgeomnav_SOURCE_DIR}, ${g4vecgeomnav_BINARY_DIR}")
    add_subdirectory(${g4vecgeomnav_SOURCE_DIR} ${g4vecgeomnav_BINARY_DIR})
  endif()
else()
  find_package(G4VecGeomNav REQUIRED)  ## 'Glue' code bridging Geant4-VecGeom navigation
endif()

set(G4VecGeomNav_LIBRARIES G4VecGeomNav::g4vecgeomnav)


#-------------------------------------------------------------------
# Locate sources for this project
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)

# define the main target :
add_executable(exampleVecGeomNav exampleVecGeomNav.cc ${sources})

target_include_directories(exampleVecGeomNav PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# we depend on G4VecGeomNav / VecGeom , optionally ROOT
target_link_libraries(exampleVecGeomNav
  PRIVATE
    ${Geant4_LIBRARIES}
    G4VecGeomNav::g4vecgeomnav
    VecGeom::vecgeom
)

# As used in G4VecGeomNav - to refine
if(WITH_ROOT)
  target_link_libraries(exampleVecGeomNav
    PRIVATE
      ROOT::Geom
      ${TGeo2VecGeom_LIBRARIES}
      ${ROOT_LIBRARIES}
  )
endif()

#----------------------------------------------------------------------------
# Copy all scripts to the build/OUTPUT directory. This is so that, after
# install, we can run the executable directly because it relies on these
# scripts being in the current working directory.
#
set(EXVECGEOMNAVIGATION_SCRIPTS
  TestNTST.gdml
  vecgeomNav.in
  run-ntst.mac
  )

foreach(_script ${EXVECGEOMNAVIGATION_SCRIPTS})
  configure_file(
    ${_script}
    ${PROJECT_BINARY_DIR}/${_script}
    COPYONLY
    )
endforeach()
