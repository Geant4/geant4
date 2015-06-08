set(_tests G4Rand*.cc G4ShootBit*.cc testRandPool.cc)
if(GEANT4_BUILD_MULTITHREADED)
	set(_tests ${_tests} testRandMT.cc testRandMarsaglia.cc testRandPerformances.cc)
endif()
GEANT4_ADD_UNIT_TESTS(${_tests} 
				 INCLUDE_DIRS global/management/include
							 global/HEPRandom/include
							${CLHEP_INCLUDE_DIRS}
				 LIBRARIES	G4global)
