# - build library
set(_TARGET _Qmaterials)
g4pytest_add_module(${_TARGET} Qmaterials.cc pyQmaterials.cc)
target_link_libraries(${_TARGET} PRIVATE G4materials)