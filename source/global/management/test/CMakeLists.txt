if(GEANT4_BUILD_MULTITHREADED)
  set(_tests testG4Cache.cc testAutoDelete.cc testG4TLSSingleton.cc)
  GEANT4_ADD_UNIT_TESTS(${_tests} 
			INCLUDE_DIRS global/management/include
			global/HEPRandom/include
			${CLHEP_INCLUDE_DIRS}
			LIBRARIES    G4global)
endif()
