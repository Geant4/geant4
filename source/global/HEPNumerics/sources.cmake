#------------------------------------------------------------------------------
# Module : G4hepnumerics
# Package: Geant4.src.G4global.G4hepnumerics
#------------------------------------------------------------------------------

#
# Define the Geant4 Module.
#
geant4_define_module(NAME G4hepnumerics
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
    G4StatAnalysis.hh
    G4StatAnalysis.icc
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
)

# List any source specific properties here
