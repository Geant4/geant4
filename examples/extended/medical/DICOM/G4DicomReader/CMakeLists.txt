#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.16...3.27)
project(G4DicomReader)

#----------------------------------------------------------------------------
# Find DCMTK installation
#
find_package(DCMTK REQUIRED)

#----------------------------------------------------------------------------
# Locate sources and headers for this project
#
file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

#----------------------------------------------------------------------------
# Add the library
#
add_library(G4DicomReader STATIC ${headers} ${sources})
target_compile_definitions(G4DicomReader PUBLIC G4_DCMTK)
target_include_directories(G4DicomReader PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>)
target_link_libraries(G4DicomReader PUBLIC ${Geant4_LIBRARIES} DCMTK::DCMTK)

add_library(G4DicomReader::G4DicomReader ALIAS G4DicomReader)


