# - build library
set(_TARGET _NISTmaterials)
g4pytest_add_module(${_TARGET} NISTmaterials.cc pyNISTmaterials.cc)
target_link_libraries(${_TARGET} PRIVATE G4materials)
