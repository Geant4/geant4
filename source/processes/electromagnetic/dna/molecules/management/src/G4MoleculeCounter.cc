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
//  Author: Mathieu Karamitros
//  Modified by Christian Velten on 10/27/2024.
//

#include "G4MoleculeCounter.hh"

#include "G4MolecularConfiguration.hh"
#include "G4Molecule.hh"
#include "G4MoleculeTable.hh"

#include <iomanip>

//------------------------------------------------------------------------------

G4MoleculeCounter::G4MoleculeCounter() : G4VUserMoleculeCounter() {}

G4MoleculeCounter::G4MoleculeCounter(G4String name)
  : G4VUserMoleculeCounter(name, MoleculeCounterType::Basic)
{}

//------------------------------------------------------------------------------

void G4MoleculeCounter::InitializeUser() {}

//------------------------------------------------------------------------------

std::unique_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
G4MoleculeCounter::BuildIndex(const G4Track* aTrack) const
{
  return std::make_unique<G4MoleculeCounterIndex>(GetMolecule(aTrack)->GetMolecularConfiguration());
}

//------------------------------------------------------------------------------

std::unique_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
G4MoleculeCounter::BuildIndex(const G4Track* aTrack, const G4StepPoint*) const
{
  return BuildIndex(aTrack);
}

//------------------------------------------------------------------------------

std::unique_ptr<G4VMoleculeCounter::G4VMoleculeCounterIndex>
G4MoleculeCounter::BuildSimpleIndex(const G4MolecularConfiguration* configuration) const
{
  return std::make_unique<G4MoleculeCounterIndex>(configuration);
}

//------------------------------------------------------------------------------
