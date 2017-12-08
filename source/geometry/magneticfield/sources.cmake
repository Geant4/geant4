#------------------------------------------------------------------------------
# sources.cmake
# Module : G4magneticfield
# Package: Geant4.src.G4geometry.G4magneticfield
#
# Sources description for a library.
# Lists the sources and headers of the code explicitely.
# Lists include paths needed.
# Lists the internal granular and global dependencies of the library.
# Source specific properties should be added at the end.
#
# Generated on : 29/9/2010
#
# $Id: sources.cmake 107113 2017-11-02 14:47:33Z gcosmo $
#
#------------------------------------------------------------------------------

# List external includes needed.
include_directories(${CLHEP_INCLUDE_DIRS})

# List internal includes needed.
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPGeometry/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/HEPRandom/include)
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4magneticfield
    HEADERS
        G4BogackiShampine23.hh
        G4BogackiShampine45.hh
        G4CachedMagneticField.hh
        G4CashKarpRKF45.hh
        G4ChargeState.hh
        G4ChordFinder.hh
        G4ChordFinder.icc
        G4ChordFinderSaf.hh
        G4ClassicalRK4.hh
        G4ConstRK4.hh
        G4DELPHIMagField.hh
        G4DoLoMcPriRK34.hh
        G4DormandPrince745.hh
        G4DormandPrinceRK56.hh
        G4DormandPrinceRK78.hh
        G4ElectricField.hh
        G4ElectroMagneticField.hh
        G4EqEMFieldWithEDM.hh
        G4EqEMFieldWithSpin.hh
        G4EqGravityField.hh
        G4EqMagElectricField.hh
        G4EquationOfMotion.hh
        G4EquationOfMotion.icc
        G4ErrorMag_UsualEqRhs.hh
        G4ExactHelixStepper.hh
        G4ExplicitEuler.hh
        G4Field.hh
        G4FieldManager.hh
        G4FieldManager.icc
        G4FieldManagerStore.hh
        G4FieldTrack.hh
        G4FieldTrack.icc
        G4FieldUtils.hh
        G4FieldUtils.icc
        G4FSALBogackiShampine45.hh
        G4FSALDormandPrince745.hh
        G4FSALIntegrationDriver.hh
        G4FSALIntegrationDriver.icc
        G4VFSALIntegrationStepper.hh
        G4VFSALIntegrationStepper.icc
        G4HarmonicPolMagField.hh
        G4HelixExplicitEuler.hh
        G4HelixHeum.hh
        G4HelixImplicitEuler.hh
        G4HelixMixedStepper.hh
        G4HelixSimpleRunge.hh
        G4ImplicitEuler.hh
        G4IntegrationDriver.hh
        G4IntegrationDriver.icc
        G4LineCurrentMagField.hh
        G4LineSection.hh
        G4MagErrorStepper.hh
        G4MagErrorStepper.icc
        G4MagHelicalStepper.hh
        G4MagHelicalStepper.icc
        G4MagIntegratorDriver.hh
        G4MagIntegratorDriver.icc
        G4MagIntegratorStepper.hh
        G4MagIntegratorStepper.icc
        G4Mag_EqRhs.hh
        G4Mag_SpinEqRhs.hh
        G4Mag_UsualEqRhs.hh
        G4MonopoleEq.hh
        G4MagneticField.hh
        G4NystromRK4.hh
        G4QuadrupoleMagField.hh
        G4RepleteEofM.hh
        G4RKG3_Stepper.hh
        G4RK547FEq1.hh
        G4RK547FEq2.hh
        G4RK547FEq3.hh
        G4SimpleHeum.hh
        G4SimpleRunge.hh
        G4TrialsCounter.hh
        G4TrialsCounter.icc
        G4TsitourasRK45.hh
        G4UniformElectricField.hh
        G4UniformGravityField.hh
        G4UniformMagField.hh
        G4VIntegrationDriver.hh
    SOURCES
        G4BogackiShampine23.cc
        G4BogackiShampine45.cc
        G4CachedMagneticField.cc
        G4CashKarpRKF45.cc
        G4ChargeState.cc
        G4ChordFinder.cc
        G4ChordFinderSaf.cc
        G4ClassicalRK4.cc
        G4ConstRK4.cc
        G4DELPHIMagField.cc
        G4DoLoMcPriRK34.cc
        G4DormandPrince745.cc
        G4DormandPrinceRK56.cc
        G4DormandPrinceRK78.cc
        G4ElectricField.cc
        G4ElectroMagneticField.cc
        G4EqEMFieldWithEDM.cc
        G4EqEMFieldWithSpin.cc
        G4EqGravityField.cc
        G4EqMagElectricField.cc
        G4EquationOfMotion.cc
        G4ErrorMag_UsualEqRhs.cc
        G4ExactHelixStepper.cc
        G4ExplicitEuler.cc
        G4Field.cc
        G4FieldManager.cc
        G4FieldManagerStore.cc
        G4FieldTrack.cc
        G4FieldUtils.cc
        G4FSALBogackiShampine45.cc
        G4FSALDormandPrince745.cc
        G4VFSALIntegrationStepper.cc
        G4HarmonicPolMagField.cc
        G4HelixExplicitEuler.cc
        G4HelixHeum.cc
        G4HelixImplicitEuler.cc
        G4HelixMixedStepper.cc
        G4HelixSimpleRunge.cc
        G4ImplicitEuler.cc
        G4LineCurrentMagField.cc
        G4LineSection.cc
        G4MagErrorStepper.cc
        G4MagHelicalStepper.cc
        G4MagIntegratorDriver.cc
        G4MagIntegratorStepper.cc
        G4Mag_EqRhs.cc
        G4Mag_SpinEqRhs.cc
        G4Mag_UsualEqRhs.cc
        G4MagneticField.cc
        G4MonopoleEq.cc
        G4NystromRK4.cc
        G4QuadrupoleMagField.cc
        G4RepleteEofM.cc
        G4RKG3_Stepper.cc
        G4RK547FEq1.cc
        G4RK547FEq2.cc
        G4RK547FEq3.cc
        G4SimpleHeum.cc
        G4SimpleRunge.cc
        G4TrialsCounter.cc
        G4TsitourasRK45.cc
        G4UniformElectricField.cc
        G4UniformGravityField.cc
        G4UniformMagField.cc
    GRANULAR_DEPENDENCIES
        G4globman
    GLOBAL_DEPENDENCIES
        G4global
    LINK_LIBRARIES
)

# List any source specific properties here

