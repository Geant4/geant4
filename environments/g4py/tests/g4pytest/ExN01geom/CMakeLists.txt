# - build library
set(_TARGET _ExN01geom)
g4pytest_add_module(${_TARGET} ExN01DetectorConstruction.cc pyExN01geom.cc)
target_link_libraries(${_TARGET} PRIVATE G4materials G4geometry G4run)
