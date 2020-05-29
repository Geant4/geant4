# - build library
set(_TARGET _ExN03geom)
g4pytest_add_module(${_TARGET} ExN03DetectorConstruction.cc ExN03DetectorMessenger.cc pyExN03geom.cc)
target_link_libraries(${_TARGET} PRIVATE G4materials G4geometry G4run)
