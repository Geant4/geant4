//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4AtomicTransition.cc,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 16 Sept 2001 Modofied according to a design iteration in the 
//              LowEnergy category
//
// -------------------------------------------------------------------

#include "G4AtomicTransition.hh"

G4AtomicTransition::G4AtomicTransition(G4int finalShell,
				       const G4std::vector<G4int>& ids,
				       const G4DataVector& energies,
				       const G4DataVector& prob)
{
  finalShellId = finalShell;
  originatingShellIds = ids;
  transitionEnergies = energies;
  transitionProbabilities = prob;
}

G4AtomicTransition::~G4AtomicTransition()
{ }

const G4std::vector<G4int>& G4AtomicTransition::OriginatingShellIds() const
{
  return  originatingShellIds;
}

const G4DataVector& G4AtomicTransition::TransitionEnergies() const
{
  return transitionEnergies;
}

const G4DataVector& G4AtomicTransition::TransitionProbabilities() const
{
  return transitionProbabilities;
}

const G4int G4AtomicTransition::FinalShellId() const
{ 
  return finalShellId;
}

G4int G4AtomicTransition::OriginatingShellId(G4int index) const
{
  return originatingShellIds[index];
}
G4double G4AtomicTransition::TransitionEnergy(G4int index) const
{
  return  transitionEnergies[index];
}
G4double G4AtomicTransition::TransitionProbability(G4int index) const
{
  return  transitionProbabilities[index];
}

