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

#ifndef DETECTORCONSTRUCTION1_HH
#define DETECTORCONSTRUCTION1_HH

#include "G4String.hh"
#include "G4Types.hh"

#include <vector>
#include <utility>

/// Define detector setup for test beam run in October 2018
/// Include the beamline elements in the simulation.
/// @param[out] aDzMap List of element names and air gap to be placed in front
/// of the element
/// @param[out] aViewpoint Targer point to be set in the visualisation
void DetectorConstruction1(
    std::vector<std::pair<G4String, G4double>> &aDzMap,
    G4double &aViewpoint);

#endif /* DETECTORCONSTRUCTION1_HH */
