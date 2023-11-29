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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
//
// File name:     G4WLSTimeGeneratorProfileDelta.cc
//
// Author:        Pedro Rodrigues, Andreia Trindade
//
// Creation date: 2006-05-07
//
// Modifications:
//
// Class Description:
//
// Class Description: End
//
// -------------------------------------------------------------------
//
//

#include "G4WLSTimeGeneratorProfileDelta.hh"
#include "Randomize.hh"
//

G4WLSTimeGeneratorProfileDelta::G4WLSTimeGeneratorProfileDelta(
  const G4String& name)
  : G4VWLSTimeGeneratorProfile(name)
{}

//

G4WLSTimeGeneratorProfileDelta::~G4WLSTimeGeneratorProfileDelta() = default;

//

G4double G4WLSTimeGeneratorProfileDelta::GenerateTime(
  const G4double time_constant)
{
  return time_constant;
}

G4double G4WLSTimeGeneratorProfileDelta::GenerateTime(
  const G4MaterialPropertiesTable*)
{
  // This method is not currently in use
  return 0.;
}
