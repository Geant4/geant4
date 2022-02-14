# - G4physicslists category build definition

geant4_global_library_target(NAME G4physicslists
  COMPONENTS
    builders/sources.cmake
	  constructors/decay/sources.cmake
	  constructors/electromagnetic/sources.cmake
	  constructors/factory/sources.cmake
	  constructors/gamma_lepto_nuclear/sources.cmake
	  constructors/hadron_elastic/sources.cmake
	  constructors/hadron_inelastic/sources.cmake
	  constructors/ions/sources.cmake
	  constructors/limiters/sources.cmake
	  constructors/stopping/sources.cmake
    lists/sources.cmake
	  util/sources.cmake)
