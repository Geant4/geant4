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

// Common implementation for histograms managers.
//
// Author: Ivana Hrivnacova, 10/08/2022  (ivana@ipno.in2p3.fr)

#ifndef G4TH2ToolsManager_h
#define G4TH2ToolsManager_h 1

#include "G4THnToolsManager.hh"

//_____________________________________________________________________________
template <>
tools::histo::p2d* G4THnToolsManager<3, tools::histo::p2d>::CreateToolsHT(
  const G4String& title,
  const std::array<G4HnDimension, 3>& bins,
  const std::array<G4HnDimensionInformation, 3>& hnInfo,
  std::string_view className);

//_____________________________________________________________________________
template <>
void G4THnToolsManager<3, tools::histo::p2d>::ConfigureToolsHT(
  tools::histo::p2d* ht,
  const std::array<G4HnDimension, 3>& bins,
  const std::array<G4HnDimensionInformation, 3>& hnInfo,
  std::string_view className);

//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<3, tools::histo::p2d>::FillHT(
  tools::histo::p2d* ht, const G4HnInformation& hnInformation, 
  std::array<G4double, 3>& value, G4double weight);


//_____________________________________________________________________________
template <>
G4bool G4THnToolsManager<3, tools::histo::p2d>::WriteOnAscii(
  std::ofstream& output);

#endif