# - build library
set(_TARGET _MedicalBeam)
g4pytest_add_module(${_TARGET} MedicalBeam.cc pyMedicalBeam.cc)
target_link_libraries(${_TARGET} PRIVATE G4run G4particles)