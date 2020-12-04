#---Adding all medical examples subdirectories explicitly 

cmake_minimum_required(VERSION 3.8...3.18)
if(${CMAKE_VERSION} VERSION_LESS 3.12)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
endif()

add_subdirectory(DICOM)
message("-- G4 Examples : The DICOM2 example requires DICOM built first.")
add_subdirectory(dna)
add_subdirectory(electronScattering)
add_subdirectory(electronScattering2)
add_subdirectory(fanoCavity)
add_subdirectory(fanoCavity2)
add_subdirectory(GammaTherapy)
