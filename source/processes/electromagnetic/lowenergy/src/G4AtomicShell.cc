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
//      Author:        Alfonso Mantero (alfonso.mantero@ge.infn.it)
// 
//      Creation date: 1/05/2001
// -------------------------------------------------------------------
#include "G4AtomicShell.hh"

// this is the constructor: U have to give some information the shell 
// wich the subshell belongs to
// the bindingenergy of the subshell u r creating, and the table of the 
// transition datas

G4AtomicShell::G4AtomicShell(G4int shellIdentifier, 
			     G4double theBindingEnergy,
			     G4std::vector<G4int> finalTransitionSubShellIdVect,
			     G4std::vector<G4double> transitionProbabilities,
			     G4std::vector<G4double> transitionEnergies) 
{
  // here u fill the variables of the subshell: the binding eneregy and 
  // the shell he belongs to

  bindingEnergy=theBindingEnergy;
  shellId=shellIdentifier;

  // these r vectors containing  datas about the transition
  // the user can see and use them thru the functions below

  finTranSubShellId = finalTransitionSubShellIdVect; 
  transProbabilities = transitionProbabilities; 
  transEnergies = transitionEnergies;
}

G4AtomicShell::G4AtomicShell()
{ }

//destructor

G4AtomicShell::~G4AtomicShell() {


  if (bindingEnergy) {
    bindingEnergy=0.0;
  }

  transProbabilities.erase(transProbabilities.begin(),transProbabilities.end());

  transEnergies.erase(transEnergies.begin(),transEnergies.end());

  finTranSubShellId.erase(finTranSubShellId.begin(),finTranSubShellId.end());


  if (shellId) {
    shellId = 0;
  }
}



G4double G4AtomicShell::BindingEnergy() const {

  return bindingEnergy;
}

const G4std::vector<G4double>& G4AtomicShell::TransitionProbabilities() const{

  return transProbabilities;
}

const G4std::vector<G4double>& G4AtomicShell::TransitionEnergies() const {

  return transEnergies;
}

const G4std::vector<G4int>& G4AtomicShell::TransSubShellIdentifiers() const{

  return finTranSubShellId;
}

G4int G4AtomicShell::ShellId() {

  return shellId;
}

G4double G4AtomicShell::TransitionEnergy(G4int index) const {

  return transEnergies[index];
}

G4int G4AtomicShell::TransitionIdentifier(G4int index) const {

  return finTranSubShellId[index];
}

G4double G4AtomicShell::TransitionProbability(G4int index) const {

  return transProbabilities[index];
}












