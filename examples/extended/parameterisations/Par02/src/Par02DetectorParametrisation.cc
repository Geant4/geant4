//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
/// \file Par02DetectorParametrisation.cc
/// \brief Implementation of the Par02DetectorParametrisation class

#include "Par02DetectorParametrisation.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02DetectorParametrisation::Par02DetectorParametrisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Par02DetectorParametrisation::~Par02DetectorParametrisation() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Par02DetectorParametrisation::GetResolution( Detector aDetector, 
                                                      Parametrisation aParam, 
                                                      G4double aMomentum ) {
  aMomentum /= GeV;  // To make sure momentum's unit is GeV
  G4double res = 1.0;
  if ( aParam == eCMS ) {
    switch ( aDetector ) {
      case Par02DetectorParametrisation::eTRACKER :
        res = 0.013;
        break;
      case Par02DetectorParametrisation::eEMCAL :
        res = std::sqrt(   std::pow( 0.03 / std::sqrt( aMomentum ), 2 )  // stochastic
                    + std::pow( 0.12 / aMomentum, 2 )          // noise
                    + std::pow( 0.003, 2 ) );                  // constant
        break;
      case Par02DetectorParametrisation::eHCAL :
        res = std::sqrt(   std::pow( 1.1 / std::sqrt( aMomentum ), 2 )   // stochastic
                    + std::pow( 0.09, 2 ) );                   // constant
        break;
    }
  } else if ( aParam == eATLAS ) {
    switch ( aDetector ) {
      case Par02DetectorParametrisation::eTRACKER :
        res = 0.01;
        break;
      case Par02DetectorParametrisation::eEMCAL :
        res = std::sqrt(   std::pow( 0.1 / std::sqrt( aMomentum ), 2 )   // stochastic
                    + std::pow( 0.0017, 2 ) );                 // constant
        break;
      case Par02DetectorParametrisation::eHCAL :
        res = std::sqrt(   std::pow( 0.55 / std::sqrt( aMomentum ), 2 )  // stochastic
                    + std::pow( 0.06, 2 ) );                   // constant
        break;
    }
  } else if ( aParam == eALEPH ) {
    switch ( aDetector ) {
      case Par02DetectorParametrisation::eTRACKER :
        res = 0.01;
        break;
      case Par02DetectorParametrisation::eEMCAL :
        res = std::sqrt(   std::pow( 0.18 / std::sqrt( aMomentum ), 2 )  // stochastic
                    + std::pow( 0.009, 2 ) );                  // constant
        break;
      case Par02DetectorParametrisation::eHCAL :
        res = 0.85 / std::sqrt( aMomentum );                   // stochastic
        break;
    }
  }
  return res;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double Par02DetectorParametrisation::GetEfficiency( Detector aDetector, 
                                                      Parametrisation /*aParam*/,
                                                      G4double /*aMomentum*/ ) {
  // For the time being, we set the efficiency to 1.0
  G4double eff = 1.0;
  switch ( aDetector ) {
    case Par02DetectorParametrisation::eTRACKER :
      eff = 1.0;
      break;
    case Par02DetectorParametrisation::eEMCAL :
      eff = 1.0;
      break;
    case Par02DetectorParametrisation::eHCAL :
      eff = 1.0;
      break;
  }
  return eff;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

