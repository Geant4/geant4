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
//      CERN, Geneva, Switzerland
//
//      File name:     G4TransitionManager.hh
//
//      Author:        Alfonso Mantero (alf@mailandnews.com)
// 
//      Creation date: 1 May 2001
//
// -------------------------------------------------------------------

#ifndef G4AtomicTransition_h
#define G4AtomicTransition_h 1

#include "G4LowEnergyUtilities.hh"
#include "G4AtomicShell.hh"
#include "globals.hh"
#include "g4std/algorithm"
#include "g4std/map"
#include "g4std/vector"
//this class is a singleton, so to be possible for only one object of this type to exist

class G4AtomicTransitionManager {

public: 

// Instance() is a function u use to create the object from outside it is the only way to do it,
// since  the creator and the destructor r protected 
  
  static G4AtomicTransitionManager* Instance();
  G4int NumberOfShells(G4int Z);
  G4double TotalRadiativeTransitionProbability(G4int Z, G4int shellId);
  G4double TotalNonRadiativeTransitionProbability(G4int Z, G4int shellId);
  G4AtomicShell* Shell(G4int Z, G4int shellIdentifier);
  // G4AtomicShell* Shell(G4int Z, G4int shellIdentifier);

protected:

  G4AtomicTransitionManager();
  ~G4AtomicTransitionManager();

private:
  
  G4LowEnergyUtilities util;
  static G4AtomicTransitionManager* instance;
  G4std::map<G4int,G4std::vector<G4AtomicShell::G4AtomicShell> > shellTable;
  void BuildTable();
};

#endif
