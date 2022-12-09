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
//
// G4PSStepChecker
#include "G4PSStepChecker.hh"

///////////////////////////////////////////////////////////////////////////////
// (Description)
//   This is a primitive scorer class for Debug.
//
// Created: 2011-03-24  Tsukasa ASO
//
///////////////////////////////////////////////////////////////////////////////

G4PSStepChecker::G4PSStepChecker(G4String name, G4int depth)
  : G4VPrimitiveScorer(name, depth)
{}

G4bool G4PSStepChecker::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  G4cout << "G4PSStepChecker:: Step identified index= " << GetIndex(aStep)
         << G4endl;
  return true;
}
