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
#ifndef DETECTORCONSTRUCTION0_HH
#define DETECTORCONSTRUCTION0_HH

#include "G4String.hh"
#include "G4Types.hh"

#include <vector>
#include <utility>

/// Add module of silicon sensor with absorbers, readout, and support, for the
/// electromagnetic calorimeter (EE).
/// @param[in] aEEid ID of module in the setup (not all modules are identical)
/// @param[out] aDzMap List of element names and air gap to be placed in front
/// of the element
/// @param[in] aAirBefore How much air should be added before the current module
/// (before the first iron)
/// @param[in] aAirMid How much air should be added in the middle (before the
/// first PCB)
void add_EE(G4int aEEid, std::vector<std::pair<G4String, G4double>> &aDzMap,
            G4double aAirBefore, G4double aAirMid);

/// Add module of scintillator sensor with absorbers, readout, and support, for
/// the hadronic calorimeter (FH).
/// @param[in] aFHid ID of module in the setup (not all modules are identical)
/// @param[out] aDzMap List of element names and air gap to be placed in front
/// of the element
/// @param[in] aAirBefore How much air should be added before the current module
/// (before the first iron)
/// @param[in] aAirMid How much air should be added in the middle (before the
/// first PCB)
void add_FH(G4int aFHid, std::vector<std::pair<G4String, G4double>> &aDzMap,
            G4double aAirBefore, G4double aAirMid);

/// Define detector setup for test beam run in October 2018
/// @param[out] aDzMap List of element names and air gap to be placed in front
/// of the element
/// @param[out] aViewpoint Targer point to be set in the visualisation
void DetectorConstruction0(
    std::vector<std::pair<G4String, G4double>> &aDzMap,
    G4double &aViewpoint);

#endif /* DETECTORCONSTRUCTION0_HH */
