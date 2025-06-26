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
// Author: Christian Velten (2025)

#include "G4MoleculeReactionCounter.hh"

//------------------------------------------------------------------------------

G4String G4MoleculeReactionCounterIndex::FormattedReactionString(const G4DNAMolecularReactionData* reactionData) const
{
  const G4MolecularConfiguration* reactant1 = reactionData->GetReactant1();
  const G4MolecularConfiguration* reactant2 = reactionData->GetReactant2();

  const std::vector<const G4MolecularConfiguration*>* products = reactionData->GetProducts();

  G4String reactionLhs = "";
  if (reactant1 != nullptr) {
    reactionLhs += reactant1->GetUserID();
    if (reactant2 != nullptr) reactionLhs += " + ";
  }
  if (reactant2 != nullptr) reactionLhs += reactant2->GetUserID();

  G4String reactionRhs = "";
  for (auto it = products->cbegin(); it != products->cend(); ++it) {
    if (*it != nullptr) {
      if (it != products->cbegin() && reactionRhs.size() > 0) reactionRhs += " + ";
      reactionRhs += (*it)->GetUserID();
    }
  }

  G4String reactionString = reactionLhs + " -> " + reactionRhs;

  return reactionString;
}

//------------------------------------------------------------------------------

G4MoleculeReactionCounter::G4MoleculeReactionCounter() : G4VUserMoleculeReactionCounter() {}

G4MoleculeReactionCounter::G4MoleculeReactionCounter(G4String name)
  : G4VUserMoleculeReactionCounter(std::move(name), MoleculeReactionCounterType::Basic)
{}

//------------------------------------------------------------------------------

void G4MoleculeReactionCounter::InitializeUser() {}

//------------------------------------------------------------------------------

std::unique_ptr<G4VMoleculeReactionCounter::G4VMoleculeReactionCounterIndex>
G4MoleculeReactionCounter::BuildSimpleIndex(const G4DNAMolecularReactionData* reactionData) const
{
  return std::make_unique<G4MoleculeReactionCounterIndex>(reactionData);
}

//------------------------------------------------------------------------------
