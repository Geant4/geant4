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
//      GEANT 4 class header file 
//
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

#include "G4DataVector.hh"
#include "globals.hh"
#include "g4std/vector"

class G4AtomicShell {

public:

 
  G4AtomicShell(G4int,G4double,const G4DataVector&,const G4DataVector&,const G4DataVector&);
 
  ~G4AtomicShell();

  //gives the binding energy of the shell whose identifier is shellId
  G4double BindingEnergy() const; 

  //the transitions considered are all the radiative transitions allowed from the 
  //shell "shellId". The identity of the final shells are contained in the 
  //vector "finalTransShellId".
  //this vector contains the probabilities of the transitions from the "shellId"
  //towards the "finalTransShellId".
  //The i-th element of the vector contains the probability of the transition
  //from the shell "shellId" towards the i-th element of the vector "finalTransShellId".
 
const G4DataVector& TransitionProbabilities() const;
  //gives the vector containing the energies of the photon emitted in the transition above 
 
const G4DataVector& TransitionEnergies() const;

  //gives a vector of the shells available from the shell "shellId" by a radiative
  //transition
const G4DataVector& TransFinalShellIdentifiers() const;
 
 //gives the starting shell 
  G4int ShellId() const;

 
  //gives the energy of the selected transition
  G4double TransitionEnergy(G4int index) const;

  //gives the identity of the final shell
  G4int TransitionIdentifier(G4int index) const;

  //gives the probability of the the selected transition
  G4double TransitionProbability(G4int index) const;

private:
  
  G4int shellId;
  G4double bindingEnergy;
G4DataVector transProbabilities;
G4DataVector transEnergies;
G4DataVector finalTranShellId;

};

#endif


