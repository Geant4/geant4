#------------------------------------------------------------------------------
# sources.cmake
# Module : G4baryons
# Package: Geant4.src.G4particles..G4baryons
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 24/9/2010
#
# $Id: sources.cmake 66892 2013-01-17 10:57:59Z gunter $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)
include_directories(${CMAKE_SOURCE_DIR}/source/materials/include)
include_directories(${CMAKE_SOURCE_DIR}/source/particles/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4baryons
    HEADERS
        G4AntiLambda.hh
	G4AntiLambdab.hh
        G4AntiLambdacPlus.hh
        G4AntiNeutron.hh
        G4AntiOmegaMinus.hh
	G4AntiOmegabMinus.hh
        G4AntiOmegacZero.hh
        G4AntiProton.hh
        G4AntiSigmaMinus.hh
        G4AntiSigmaPlus.hh
        G4AntiSigmaZero.hh
        G4AntiSigmacPlus.hh
        G4AntiSigmacPlusPlus.hh
        G4AntiSigmacZero.hh
	G4AntiSigmabMinus.hh
	G4AntiSigmabPlus.hh
	G4AntiSigmabZero.hh
        G4AntiXiMinus.hh
        G4AntiXiZero.hh
        G4AntiXicPlus.hh
        G4AntiXicZero.hh
	G4AntiXibMinus.hh 
	G4AntiXibZero.hh
        G4BaryonConstructor.hh
        G4Lambda.hh
	G4Lambdab.hh
        G4LambdacPlus.hh
        G4Neutron.hh
        G4OmegaMinus.hh
        G4OmegabMinus.hh
        G4OmegacZero.hh
        G4Proton.hh
        G4SigmaMinus.hh
        G4SigmaPlus.hh
        G4SigmaZero.hh
        G4SigmacPlus.hh
        G4SigmacPlusPlus.hh
        G4SigmacZero.hh
	G4SigmabMinus.hh
	G4SigmabPlus.hh
	G4SigmabZero.hh
        G4XiMinus.hh
        G4XiZero.hh
        G4XicPlus.hh
        G4XicZero.hh
	G4XibMinus.hh 
	G4XibZero.hh
    SOURCES
        G4AntiLambda.cc
	G4AntiLambdab.cc
        G4AntiLambdacPlus.cc
        G4AntiNeutron.cc
        G4AntiOmegaMinus.cc
        G4AntiOmegabMinus.cc
        G4AntiOmegacZero.cc
        G4AntiProton.cc
        G4AntiSigmaMinus.cc
        G4AntiSigmaPlus.cc
        G4AntiSigmaZero.cc
        G4AntiSigmacPlus.cc
        G4AntiSigmacPlusPlus.cc
        G4AntiSigmacZero.cc
	G4AntiSigmabMinus.cc
	G4AntiSigmabPlus.cc
	G4AntiSigmabZero.cc
        G4AntiXiMinus.cc
        G4AntiXiZero.cc
        G4AntiXicPlus.cc
        G4AntiXicZero.cc
	G4AntiXibMinus.cc
	G4AntiXibZero.cc
        G4BaryonConstructor.cc
        G4Lambda.cc
	G4Lambdab.cc
        G4LambdacPlus.cc
        G4Neutron.cc
        G4OmegaMinus.cc
        G4OmegabMinus.cc
        G4OmegacZero.cc
        G4Proton.cc
        G4SigmaMinus.cc
        G4SigmaPlus.cc
        G4SigmaZero.cc
        G4SigmacPlus.cc
        G4SigmacPlusPlus.cc
        G4SigmacZero.cc
	G4SigmabMinus.cc
	G4SigmabPlus.cc
	G4SigmabZero.cc
        G4XiMinus.cc
        G4XiZero.cc
        G4XicPlus.cc
        G4XicZero.cc
	G4XibMinus.cc
	G4XibZero.cc
    GRANULAR_DEPENDENCIES
        G4globman
        G4materials
        G4partman
    GLOBAL_DEPENDENCIES
        G4global
        G4materials
    LINK_LIBRARIES
)

# List any source specific properties here

