#------------------------------------------------------------------------------
# sources.cmake
# Module : G4hepnumerics
# Package: Geant4.src.G4global.G4hepnumerics
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
include_directories(${CMAKE_SOURCE_DIR}/source/global/management/include)

#
# Define the Geant4 Module.
#
include(Geant4MacroDefineModule)
GEANT4_DEFINE_MODULE(NAME G4hepnumerics
    HEADERS
        G4AnalyticalPolSolver.hh
        G4ChebyshevApproximation.hh
        G4ConvergenceTester.hh
        G4DataInterpolation.hh
        G4GaussChebyshevQ.hh
        G4GaussHermiteQ.hh
        G4GaussJacobiQ.hh
        G4GaussLaguerreQ.hh
        G4GaussLegendreQ.hh
        G4Integrator.hh
        G4Integrator.icc
        G4JTPolynomialSolver.hh
        G4PolynomialSolver.hh
        G4PolynomialSolver.icc
        G4SimpleIntegration.hh
        G4SimplexDownhill.hh
        G4SimplexDownhill.icc
        G4StatDouble.hh
        G4VGaussianQuadrature.hh
    SOURCES
        G4AnalyticalPolSolver.cc
        G4ChebyshevApproximation.cc
        G4ConvergenceTester.cc
        G4DataInterpolation.cc
        G4GaussChebyshevQ.cc
        G4GaussHermiteQ.cc
        G4GaussJacobiQ.cc
        G4GaussLaguerreQ.cc
        G4GaussLegendreQ.cc
        G4JTPolynomialSolver.cc
        G4SimpleIntegration.cc
        G4StatDouble.cc
        G4VGaussianQuadrature.cc
    GRANULAR_DEPENDENCIES
        G4globman
    GLOBAL_DEPENDENCIES
    LINK_LIBRARIES
)

# List any source specific properties here

