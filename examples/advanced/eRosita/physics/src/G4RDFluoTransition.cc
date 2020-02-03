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
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 16 Sept 2001  EG  Modified according to a design iteration in the 
//                   LowEnergy category
//
// -------------------------------------------------------------------

#include "G4RDFluoTransition.hh"

G4RDFluoTransition::G4RDFluoTransition(G4int finalShell,
				       const std::vector<G4int>& ids,
				       const G4DataVector& energies,
				       const G4DataVector& prob)
  :finalShellId(finalShell),
   originatingShellIds(ids),
   transitionEnergies(energies),
   transitionProbabilities(prob)
{ }

G4RDFluoTransition::~G4RDFluoTransition()
{ }

const std::vector<G4int>& G4RDFluoTransition::OriginatingShellIds() const
{
  return  originatingShellIds;
}

const G4DataVector& G4RDFluoTransition::TransitionEnergies() const
{
  return transitionEnergies;
}

const G4DataVector& G4RDFluoTransition::TransitionProbabilities() const
{
  return transitionProbabilities;
}

G4int G4RDFluoTransition::FinalShellId() const
{ 
  return finalShellId;
}

G4int G4RDFluoTransition::OriginatingShellId(G4int index) const
{
  return originatingShellIds[index];
}
G4double G4RDFluoTransition::TransitionEnergy(G4int index) const
{
  return  transitionEnergies[index];
}
G4double G4RDFluoTransition::TransitionProbability(G4int index) const
{
  return  transitionProbabilities[index];
}

