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
#include "DetectorConstruction1.hh"
#include "DetectorConstruction0.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction1(
    std::vector<std::pair<G4String, G4double>> &aDzMap,
    G4double &aViewpoint) {

  // map means: position this material starting at z, where:
  // z += <second>
  // z = z0+ 0.5 * thickness[<first>]

  // first define beamline elements

  // WChambUpstream
  // world (and beam) starts at -45 m
  aDzMap.push_back(std::make_pair("DWC", 12.2825 * CLHEP::m)); // at -32.7175 m
  aDzMap.push_back(std::make_pair("CK3", 0.0 * CLHEP::m));     // at -30.69 m
  aDzMap.push_back(std::make_pair("DWC", 56.2 * CLHEP::cm));   // at -30.073 cm
  aDzMap.push_back(std::make_pair("DWC", 101.5 * CLHEP::cm));  // at -2900.3 cm
  aDzMap.push_back(std::make_pair("Scintillator_thin", 23.9 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("Scintillator_thin", 13.9 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("Scintillator_thin", 13.9 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("Scintillator_thin", 103.8 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("DWC", 34.8 * CLHEP::cm)); // at -2700.5 cm
  // HaloCounters not implemented (as offset in X or Y)

  // WChambDown
  aDzMap.push_back(std::make_pair("DWC", 18.15 * CLHEP::m)); // at -880 cm
  aDzMap.push_back(std::make_pair("DWC", 7.145 * CLHEP::m)); // at -160 cm

  // S5
  aDzMap.push_back(std::make_pair("Scintillator", 1.354 * CLHEP::m));
  // S6
  aDzMap.push_back(std::make_pair("Scintillator", 14 * CLHEP::cm));

  // next use already defined HGCal configuration

  DetectorConstruction0(aDzMap, aViewpoint);
}
