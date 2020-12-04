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
#include "DetectorConstruction0.hh"
#include "CLHEP/Units/SystemOfUnits.h"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void add_EE(G4int aEEid, std::vector<std::pair<G4String, G4double>> &aDzMap,
            G4double aAirBefore, G4double aAirMid) {
  aDzMap.push_back(std::make_pair("Fe_absorber_EE",
                                  aAirBefore)); // aAirBefore AIR + 0.3 mm Fe
  aDzMap.push_back(std::make_pair("Pb_absorber_EE", 0)); // 0 + 4.9mm Pb
  aDzMap.push_back(std::make_pair("Fe_absorber_EE", 0)); // 0 + 0.3 mm Fe
  aDzMap.push_back(std::make_pair("PCB", aAirMid));      // aAirMid AIR +1.3mm
  aDzMap.push_back(std::make_pair("Cu_baseplate_175um", 0.)); // 0 + 0.175mm Cu
  aDzMap.push_back(std::make_pair("Si_wafer", 0.));           // 0.3 mm
  aDzMap.push_back(std::make_pair("Cu_baseplate_25um", 0.));  //  0 + 0.025mm Cu
  aDzMap.push_back(std::make_pair("Kapton_layer", 0.));       // 0.075 mm
  if (aEEid == 11 || aEEid == 12)
    aDzMap.push_back(std::make_pair("Cu_baseplate", 0.)); // 1.2 mm
  aDzMap.push_back(std::make_pair("CuW_baseplate", 0.));  // 1.2 mm
  if (aEEid == 13)
    aDzMap.push_back(std::make_pair("CuW_baseplate_550um", 0.)); // 0.55 mm
  if (aEEid == 14)
    aDzMap.push_back(std::make_pair("CuW_baseplate_610um", 0.)); // 0.61 mm
  aDzMap.push_back(std::make_pair("Cu_absorber_EE", 0.));        // 6 mm
  if (aEEid == 13)
    aDzMap.push_back(std::make_pair("CuW_baseplate_610um", 0.)); // 0.61 mm
  if (aEEid == 14)
    aDzMap.push_back(std::make_pair("CuW_baseplate_710um", 0.)); // 0.71 mm
  aDzMap.push_back(std::make_pair("CuW_baseplate", 0.));         // 1.2 mm
  if (aEEid == 11 || aEEid == 12)
    aDzMap.push_back(std::make_pair("Cu_baseplate", 0.));     // 1.2 mm
  aDzMap.push_back(std::make_pair("Kapton_layer", 0.));       // 0.075 mm
  aDzMap.push_back(std::make_pair("Cu_baseplate_25um", 0.));  //  0 + 0.025mm Cu
  aDzMap.push_back(std::make_pair("Si_wafer", 0.));           // 0.3 mm
  aDzMap.push_back(std::make_pair("Cu_baseplate_175um", 0.)); // 0 + 0.175mm Cu
  aDzMap.push_back(std::make_pair("PCB", 0));                 // 1.3 mm
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void add_FH(G4int aFHid, std::vector<std::pair<G4String, G4double>> &aDzMap,
            G4double aAirBefore, G4double aAirMid) {
  std::string layout = "";
  if (aFHid < 10)
    layout = "_DAISY";
  if (!(aFHid == 1 || aFHid == 7))
    aDzMap.push_back(std::make_pair("Fe_absorber_FH",
                                    aAirBefore)); // aAirBefore AIR + 40 mm Fe
  aDzMap.push_back(
      std::make_pair("PCB" + layout, aAirMid)); // aAirMid AIR + 1.3 mm
  aDzMap.push_back(
      std::make_pair("Cu_baseplate_175um" + layout, 0.));    // 0.175mm Cu
  aDzMap.push_back(std::make_pair("Si_wafer" + layout, 0.)); // 0.3 mm
  if (aFHid == 5) {
    aDzMap.push_back(std::make_pair("PCB_thin" + layout, 0)); // 1.2 mm
  }
  if (aFHid != 5) {
    aDzMap.push_back(
        std::make_pair("Cu_baseplate_25um" + layout, 0.));         // 0.025mm Cu
    aDzMap.push_back(std::make_pair("Kapton_layer" + layout, 0.)); // 0.075 mm
  }
  if (aFHid == 6) {
    aDzMap.push_back(
        std::make_pair("Cu_baseplate_25um" + layout, 0.));         // 0.025mm Cu
    aDzMap.push_back(std::make_pair("Kapton_layer" + layout, 0.)); // 0.075 mm
  }
  if (aFHid != 5 && aFHid < 9)
    aDzMap.push_back(std::make_pair("Cu_baseplate" + layout, 0.)); // 1.2 mm
  if (aFHid == 9 || aFHid == 10)
    aDzMap.push_back(std::make_pair("CuW_baseplate" + layout, 0.)); // 1.2 mm
  if (aFHid != 10)
    aDzMap.push_back(std::make_pair("Cu_baseplate" + layout, 0.)); // 1.2 mm
  aDzMap.push_back(std::make_pair("Cu_absorber_FH", 0.));          // 6 * mm
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction0(
    std::vector<std::pair<G4String, G4double>> &aDzMap,
    G4double &aViewpoint) {
  aViewpoint = 0.75 * CLHEP::m;

  // map means: position this material starting at z, where:
  // z += <second>
  // z = z0+ 0.5 * thickness[<first>]

  // world (and beam) starts at -45 m

  G4double firstOffset = 2.6 * CLHEP::cm;
  // if no beamline is present - shift by beamline length 45.0015 m
  if (aDzMap.size() == 0) {
    firstOffset += 45.0015 * CLHEP::m;
  }

  aDzMap.push_back(std::make_pair("Al_case_thick", firstOffset)); // 5 mm Al
  aDzMap.push_back(std::make_pair("Al_case", 0));              // 2.1 mm Al

  // EE1
  add_EE(1, aDzMap, 119.7 * CLHEP::mm, 4.7 * CLHEP::mm);

  // EE2
  add_EE(2, aDzMap, 7.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE3
  add_EE(3, aDzMap, 7.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE4
  add_EE(4, aDzMap, 8.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE5
  add_EE(5, aDzMap, 8.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE6
  add_EE(6, aDzMap, 8.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE7
  add_EE(7, aDzMap, 6.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE8
  add_EE(8, aDzMap, 6.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE9
  add_EE(9, aDzMap, 6.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE10
  add_EE(10, aDzMap, 6.7 * CLHEP::mm, 3.7 * CLHEP::mm);

  // EE11
  add_EE(11, aDzMap, 6.7 * CLHEP::mm, 5.5 * CLHEP::mm);

  // EE12
  add_EE(12, aDzMap, 9.5 * CLHEP::mm, 5.5 * CLHEP::mm);

  // EE13
  add_EE(13, aDzMap, 9.5 * CLHEP::mm, 3.145 * CLHEP::mm);

  // EE14
  add_EE(14, aDzMap, 10.09 * CLHEP::mm, 3.095 * CLHEP::mm);

  aDzMap.push_back(std::make_pair("Steel_case_thick", 0)); // 40 mm
  aDzMap.push_back(std::make_pair("Al_case", 44 * CLHEP::mm));    // 2.1 mm

  // beginning of FH
  aDzMap.push_back(std::make_pair("Steel_case", 0)); // 9 mm

  // FH1
  add_FH(1, aDzMap, 0, 8.8 * CLHEP::mm);

  // FH2
  add_FH(2, aDzMap, 8 * CLHEP::mm, 8.8 * CLHEP::mm);

  // FH3
  add_FH(3, aDzMap, 3 * CLHEP::mm, 13.8 * CLHEP::mm);

  // FH4
  add_FH(4, aDzMap, 5 * CLHEP::mm, 12.8 * CLHEP::mm);

  // FH5
  add_FH(5, aDzMap, 7 * CLHEP::mm, 9.8 * CLHEP::mm);

  // FH6
  add_FH(6, aDzMap, 6 * CLHEP::mm, 10.7 * CLHEP::mm);

  // cases
  aDzMap.push_back(std::make_pair("Steel_case", 4 * CLHEP::mm));      // 9 mm
  aDzMap.push_back(std::make_pair("Fe_absorber_FH", 36 * CLHEP::mm)); // 40 mm
  aDzMap.push_back(std::make_pair("Steel_case", 52 * CLHEP::mm));     // 9 mm
  // FH7
  add_FH(7, aDzMap, 0, 8.8 * CLHEP::mm);

  // FH8
  add_FH(8, aDzMap, 7 * CLHEP::mm, 16.8 * CLHEP::mm);

  // FH9
  add_FH(9, aDzMap, 9 * CLHEP::mm, 14.8 * CLHEP::mm);

  // FH10
  add_FH(10, aDzMap, 10 * CLHEP::mm, 18 * CLHEP::mm);

  // FH11
  add_FH(11, aDzMap, 8 * CLHEP::mm, 17 * CLHEP::mm);

  // FH12
  add_FH(12, aDzMap, 7 * CLHEP::mm, 17 * CLHEP::mm);

  aDzMap.push_back(std::make_pair("Steel_case", 29 * CLHEP::mm));

  // AHCAL
  aDzMap.push_back(std::make_pair("Fe_absorber_AHCAL", 50.0 * CLHEP::cm));
  for (int l = 0; l < 39; l++) {
    aDzMap.push_back(std::make_pair("Al_absorber_AHCAL", 0.5 * CLHEP::cm));
    aDzMap.push_back(std::make_pair("AHCAL_SiPM_2x2HUB", 0.));
    aDzMap.push_back(std::make_pair("Al_absorber_AHCAL", 0.));
    aDzMap.push_back(std::make_pair("Fe_absorber_AHCAL", 0.5 * CLHEP::cm));
  }
  aDzMap.push_back(std::make_pair("Fe_absorber_AHCAL", 1.1 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("Fe_absorber_AHCAL", 1.1 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("Al_absorber_AHCAL", 0.5 * CLHEP::cm));
  aDzMap.push_back(std::make_pair("AHCAL_SiPM_2x2HUB", 0.));
  aDzMap.push_back(std::make_pair("Al_absorber_AHCAL", 0.));
  aDzMap.push_back(std::make_pair("Fe_absorber_AHCAL", 0.5 * CLHEP::cm));
}
