#----------------------------------------------------------------------------
# Build G4DicomCore library for shared functionality between the two examples
#----------------------------------------------------------------------------

# Locate sources and headers for this project
# NB: headers are included so they will show up in IDEs
#
file(GLOB sources ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${CMAKE_CURRENT_SOURCE_DIR}/include/*.hh)

add_library(G4DicomCore STATIC ${headers} ${sources})
target_compile_definitions(G4DicomCore PRIVATE $<$<BOOL:${G4DICOM_USE_HEAD}>:DICOM_USE_HEAD>)
target_include_directories(G4DicomCore PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_link_libraries(G4DicomCore 
  PUBLIC
    ${Geant4_LIBRARIES}
    $<$<TARGET_EXISTS:G4DicomReader::G4DicomReader>:G4DicomReader::G4DicomReader>)

add_library(G4DicomCore::G4DicomCore ALIAS G4DicomCore)
