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
// $Id: G4FluoTransition.cc,v 1.2 ????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
// 16 Sept 2001  EG  Modified according to a design iteration in the 
//                   LowEnergy category
//
// -------------------------------------------------------------------

#include "G4FluoTransition.hh"

G4FluoTransition::G4FluoTransition(G4int finalShell,
				       const G4std::vector<G4int>& ids,
				       const G4DataVector& energies,
				       const G4DataVector& prob)
  :finalShellId(finalShell),
   originatingShellIds(ids),
   transitionEnergies(energies),
   transitionProbabilities(prob)
{ }

G4FluoTransition::~G4FluoTransition()
{ }

const G4std::vector<G4int>& G4FluoTransition::OriginatingShellIds() const
{
  return  originatingShellIds;
}

const G4DataVector& G4FluoTransition::TransitionEnergies() const
{
  return transitionEnergies;
}

const G4DataVector& G4FluoTransition::TransitionProbabilities() const
{
  return transitionProbabilities;
}

const G4int G4FluoTransition::FinalShellId() const
{ 
  return finalShellId;
}

G4int G4FluoTransition::OriginatingShellId(G4int index) const
{
  return originatingShellIds[index];
}
G4double G4FluoTransition::TransitionEnergy(G4int index) const
{
  return  transitionEnergies[index];
}
G4double G4FluoTransition::TransitionProbability(G4int index) const
{
  return  transitionProbabilities[index];
}

