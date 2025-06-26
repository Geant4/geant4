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
//  G4VMoleculeCounter.cc
//  Geant4
//
//  Created by Mathieu Karamitros on 02/11/2016.
//  Modified by Christian Velten on 10/27/2024.
//
//

#include "G4VMoleculeCounter.hh"

#include "G4MoleculeCounterTemplates.hh"
#include "G4ios.hh"

//------------------------------------------------------------------------------

G4VMoleculeCounter::G4VMoleculeCounter() : G4VMoleculeCounterInternalBase() {}

G4VMoleculeCounter::G4VMoleculeCounter(const G4String& name, MoleculeCounterType type)
  : G4VMoleculeCounterInternalBase(name), fType(type)
{}

//------------------------------------------------------------------------------

void G4VMoleculeCounter::SetSensitiveToStepping(G4bool flag)
{
  if (fType == MoleculeCounterType::Basic && flag) {
    G4ExceptionDescription errMsg;
    errMsg << "Cannot set a molecule counter of type 'Basic' to be sensitive to stepping!"
           << G4endl;
    G4Exception("G4VMoleculeCounter::SetSensitiveToStepping", "NOT_ALLOWED", FatalException, errMsg);
  }
  fSensitiveToStepping = flag;
}

//------------------------------------------------------------------------------

G4bool G4VMoleculeCounter::IsReactantIgnored(const G4MoleculeDefinition* molecule) const
{
  return G4::MoleculeCounter::Contains(fIgnoredMolecules, molecule);
}

//------------------------------------------------------------------------------

G4bool G4VMoleculeCounter::IsReactantIgnored(const G4MolecularConfiguration* reactant) const
{
  return G4::MoleculeCounter::Contains(fIgnoredReactants, reactant);
}
