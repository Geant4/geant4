# - build library
set(_TARGET _EMSTDpl)
g4pytest_add_module(${_TARGET} PhysicsListEMstd.cc pyEMSTDpl.cc)
target_link_libraries(${_TARGET} PRIVATE G4particles G4processes G4run)
