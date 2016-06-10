#------------------------------------------------------------------------------
# sources.cmake
# Module : 
# Package: Geant4.src.G4global.
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 89393 2015-04-09 07:45:40Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#

# The following are needed only in
# multi-threaded builds
set (_mtheaders "")
set (_mtsrcs "")
if(GEANT4_BUILD_MULTITHREADED)
  set(_mtheaders ${_mtheaders} 
	G4MTHepRandom.hh 
	G4MTHepRandom.icc
	G4MTRandBit.hh
	G4MTRandBit.icc
	G4MTRandExponential.hh
	G4MTRandExponential.icc
	G4MTRandFlat.hh
	G4MTRandFlat.icc
	G4MTRandGamma.hh
	G4MTRandGamma.icc
	G4MTRandGauss.hh
	G4MTRandGauss.icc
	G4MTRandGaussQ.hh
	G4MTRandGaussQ.icc
	G4MTRandGeneral.hh
	G4MTRandGeneral.icc )
  set(_mtsrcs ${_mtsrcs}
	G4MTHepRandom.cc 
	G4MTRandBit.cc 
	G4MTRandExponential.cc
	G4MTRandFlat.cc
	G4MTRandGamma.cc
	G4MTRandGauss.cc
	G4MTRandGaussQ.cc
	G4MTRandGeneral.cc )
endif()

include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4heprandom 
    HEADERS
        G4Poisson.hh
        G4RandomDirection.hh
        G4RandomTools.hh
        Randomize.hh
	G4UniformRandPool.hh
	${_mtheaders}
    SOURCES
	${_mtsrcs}
        G4Poisson.cc
	G4UniformRandPool.cc
    GRANULAR_DEPENDENCIES
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

