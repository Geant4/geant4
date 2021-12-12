# - build library
set(_TARGET _Qgeom)
g4pytest_add_module(${_TARGET} QDetectorConstruction.cc pyQgeom.cc)
target_link_libraries(${_TARGET} PRIVATE G4materials G4geometry G4run)
