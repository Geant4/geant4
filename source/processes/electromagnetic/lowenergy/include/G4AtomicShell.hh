// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// -------------------------------------------------------------------
//      GEANT 4 class header file 
//
//      For information related to this code contact:
//      Geant4 Collaboration
//
//      File name:     G4AtomicShell.hh
//
//      Author:        Alfonso Mantero (alfonso.mantero@ge.infn.it)
// 
//      Creation date: 1 May 2001
//
// -------------------------------------------------------------------

#ifndef G4AtomicShell_h 
#define G4AtomicShell_h 1


#include "globals.hh"
#include "g4std/vector"

class G4AtomicShell {

public:

  G4AtomicShell(G4int shellIdentifier, 
		G4double theBindingEnergy,
		G4std::vector<G4int> finalTransitionSubShellIdVect,
		G4std::vector<G4double> transitionProbabilities,
		G4std::vector<G4double> transitionEnergies);

  G4AtomicShell();

  // destructor
  ~G4AtomicShell();

  //these functions return shell datas 

  //gives the binding energy of the subshell
  G4double BindingEnergy() const; 

  //gives a vector of transition probabilities
  const G4std::vector<G4double>& TransitionProbabilities() const; 

  //gives a vector of transition energies
  const G4std::vector<G4double>& TransitionEnergies() const;

  //gives a vector of the final subshells of the transitions
  const G4std::vector<G4int>& TransSubShellIdentifiers() const;

  //gives the shell which the subshell belongs to
  G4int ShellId();

  //gives the energy if the selected transition
  G4double TransitionEnergy(G4int index) const;

  //gives the final subshell ID
  G4int TransitionIdentifier(G4int index) const;

  //gives the probability of the the selected transition
  G4double TransitionProbability(G4int index) const;

private:
  // queste sono le variabili (vettori) dove vengono messi i dati
  G4int shellId;
  G4double bindingEnergy;
  G4std::vector<G4double> transProbabilities;
  G4std::vector<G4double> transEnergies;
  G4std::vector<G4int> finTranSubShellId;

};

#endif


