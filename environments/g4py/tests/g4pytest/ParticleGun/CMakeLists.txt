# - build library
set(_TARGET _ParticleGun)
g4pytest_add_module(${_TARGET} ParticleGunAction.cc pyParticleGun.cc)
target_link_libraries(${_TARGET} PRIVATE G4run G4particles)
