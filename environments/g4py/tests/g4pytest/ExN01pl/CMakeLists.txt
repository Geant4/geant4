# - build library
set(_TARGET _ExN01pl)
g4pytest_add_module(${_TARGET} ExN01PhysicsList.cc pyExN01pl.cc)
target_link_libraries(${_TARGET} PRIVATE G4particles G4processes G4run)
