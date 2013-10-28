if(__HEADERS_AND_SOURCES_INCLUDED)
  return()
endif()
set(__HEADERS_AND_SOURCES_INCLUDED TRUE)


function(HEADERS_AND_SOURCES _headers _sources)
#
#----------------------------------------------------------------------------
# Find Geant4 package, no UI and Vis drivers activated
#
#find_package(Geant4 REQUIRED)

#----------------------------------------------------------------------------
# setup include directories
#
include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include
                    ${Geant4_INCLUDE_DIR})

#----------------------------------------------------------------------------
#
# sources and headers for this project
#
file(GLOB prj_sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB prj_headers ${PROJECT_SOURCE_DIR}/include/*.hh)

message("Prj Dir: ${PROJECT_SOURCE_DIR}")
message("prj_sources: ${prj_sources}")

#
#  sources from shared code
#
file(GLOB g4app_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/src/*.cc)
file(GLOB g4app_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/g4app/include/*.hh)

# find additional ROOT package
#
#find_package(ROOT QUIET)
#if(NOT ROOT_FOUND)
#  message(STATUS "G4 TESTS: ROOT package not found. ")  
#  return()
#endif()

#
# additional headers and spources from ROOT (if any)
#
include_directories(${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/include
		    ${ROOTSYS}/include)
file(GLOB rootana_sources 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootanalysis/src/*.cc)
file(GLOB rootana_headers 
  ${PROJECT_SOURCE_DIR}/../test23/shared/rootana/include/*.hh)

#
# finally, put all headers and sources together 
#
list(APPEND sources ${prj_sources} ${g4app_sources} ${rootana_sources})
list(APPEND headers ${prj_headers} ${g4app_headers} ${rootana_headers})

set(${_sources} ${sources})
set(${_headers} ${headers})

endfunction(HEADERS_AND_SOURCES _headers _sources)
