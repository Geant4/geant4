// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
//      For information related to this code contact:
//      Geant4 Collaboration
//
//      File name:     G4AtomicShell.cc
//
//      Authors:        Alfonso Mantero (alfonso.mantero@ge.infn.it)
//                      Elena Guardincerri (elena.guardincerri@ge.infn.it) 
//
//      Creation date: 1/05/2001
// -------------------------------------------------------------------
#include "G4AtomicShell.hh"

G4AtomicShell::G4AtomicShell(G4int id, G4double energy, 
			     const G4DataVector& prob,
			     const G4DataVector& energies,
			     const G4DataVector& shells)
{
  shellId = id;
  bindingEnergy = energy;
  transProbabilities = prob;
  transEnergies = energies;
  finalTranShellId = shells;

} 

G4AtomicShell::~G4AtomicShell()
{;}

G4double G4AtomicShell::BindingEnergy() const {

  return bindingEnergy;
}

const G4DataVector& G4AtomicShell::TransitionProbabilities() const{
 
  return transProbabilities;
}

const G4DataVector& G4AtomicShell::TransitionEnergies() const {

  return transEnergies;
}

const G4DataVector& G4AtomicShell::TransFinalShellIdentifiers() const{

  return finalTranShellId;
}

G4int G4AtomicShell::ShellId() const{

  return shellId;
}

G4double G4AtomicShell::TransitionEnergy(G4int index) const {

  return transEnergies[index];
}

G4int G4AtomicShell::TransitionIdentifier(G4int index) const {

  return (G4int)finalTranShellId[index];
}

G4double G4AtomicShell::TransitionProbability(G4int index) const {

  return transProbabilities[index];
}





