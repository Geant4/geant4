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
// -------------------------------------------------------------------
//      GEANT 4 class implementation file 
//
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

//const G4std::vector<G4double>& G4AtomicShell::TransitionProbabilities() const{
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





