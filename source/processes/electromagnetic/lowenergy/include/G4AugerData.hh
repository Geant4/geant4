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
// $Id: G4AugerData.hh
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Alfonso Mantero (Alfonso.Mantero@ge.infn.it)
//
// History:
// -----------
//  2 June 2002 First committed to cvs
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// This Class loads and stores all the information of auger effect (shellIds, 
// probabilities and  energies of the electrons emitted) 
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4AUGERDATA_HH
#define G4AUGERDATA_HH 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"
#include "G4AugerTransition.hh"

class G4DataVector;

class G4AugerData
{
public:

  G4AugerData();

  ~G4AugerData();

  // The method returns the number of shells in wich a 
  // vacancy can be filled by a NON-radiative transition, given the atomic number
  size_t NumberOfVacancies(G4int Z) const;

  // Given the index of the vacancy (and the atomic number Z) returns its identity
  G4int VacancyId(G4int Z, G4int vacancyIndex) const;
  
  // Given the index of a vacancy in the atom with the atomc number Z, returns the number of
  //shells starting from wich an electron can fill the vacancy
  size_t NumberOfTransitions(G4int Z, G4int vacancyIndex) const;

  // Given the atomic number Z, the Index of the initial vacancy shell 
  // and the index of the starting shell for the 
  // transition, returns the identity of the shell originating the electron transition
  G4int StartShellId(G4int Z, G4int initialVacancyIndex, G4int transitionShellIndex) const;

  // Given the atomic number , the indexes of the starting, the auger originating shell, 
  // and the transition shell Id, returns the transition energy
  G4double StartShellEnergy(G4int Z, G4int vacancyIndex, G4int transitionId, G4int augerIndex) const;

  // Given the atomic number, the  index of the starting shell, the auger originating shells, 
  // and the transition shell Id, returns the transition probability
  G4double StartShellProb(G4int Z, G4int vacancyIndex,G4int transitionId,G4int augerIndex) const;

  // Given the atomic number, the index of the starting vacancy shell and the transition shell Id,
  // returns the number of shells wich an auger electron can come from.
  size_t NumberOfAuger(G4int Z, G4int initIndex, G4int vacancyId) const;

  // Given the atomic number, th index of the starting and the auger originating shell, 
  // and the transition shell Id, returns the ager originating shell Id
  size_t AugerShellId(G4int Z, G4int vacancyIndex, G4int transId, G4int augerIndex) const;

  G4std::vector<G4AugerTransition> LoadData(G4int Z);

  void BuildAugerTransitionTable();

  void PrintData(G4int Z);



  // Given the atomic number and the vacancy intial shell index  returns 
  // the AugerTransition object related to that shell

  G4AugerTransition* GetAugerTransition(G4int Z, G4int vacancyShellIndex);
  
  // Given the atomic number returns a vector of possible AugerTransition objects
  G4std::vector<G4AugerTransition>* GetAugerTransitions(G4int Z);

private:

  // G4std::map<G4int,G4DataVector*,G4std::less<G4int> > idMap;

  typedef G4std::map<G4int,G4std::vector<G4AugerTransition>,G4std::less<G4int> > trans_Table;
   trans_Table augerTransitionTable;

  /*
  G4std::map<G4int,G4std::map<G4Int,G4DataVector*,G4std::less<G4int> >,G4std::less<G4int> > transProbabilityMap;
  G4std::map<G4int,G4std::map<G4Int,G4DataVector*,G4std::less<G4int> >,G4std::less<G4int> > transAugerIdMap;
  */

  G4std::vector<G4int> nInitShells;
  G4std::vector<G4int> numberOfVacancies;
  
};

#endif





