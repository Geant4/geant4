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
// $Id: G4FluoData.hh,v 1.1 ?????
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Author: Elena Guardincerri (Elena.Guardincerri@ge.infn.it)
//
// History:
// -----------
//  16 Sept 2001 First committed to cvs
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Fluorescence data set: shell identifiers, transition probabilities, 
// transition energies
// Further documentation available from http://www.ge.infn.it/geant4/lowE

// -------------------------------------------------------------------

#ifndef G4FLUODATA_HH
#define G4FLUODATA_HH 1

#include "globals.hh"
#include "g4std/vector"
#include "g4std/map"


class G4DataVector;

class G4FluoData
{
public:

  G4FluoData();

  ~G4FluoData();

  // The method returns the number of shells in wich a 
  // vacancy can be filled by a radiative transition
  size_t NumberOfVacancies() const;

  // Given the index of the vacancy returns its identity
  G4int VacancyId(G4int vacancyIndex) const;
  
  // Given the index of a vacancy returns the number of
  //shells starting from wich an electrons can fill the vacancy
  size_t NumberOfTransitions(G4int vacancyIndex) const;

  // Given the indexes of the starting and final shells for the 
  // transition, returns the identity of the starting one
  G4int StartShellId(G4int initIndex,G4int vacancyIndex);

  // Given the indexes of the starting and final shells for the 
  // transition, returns the transition energy
  G4double StartShellEnergy(G4int initIndex,G4int vacancyIndex);

  // Given the indexes of the starting and final shells for the 
  // transition, returns the probability of this transition
  G4double StartShellProb(G4int initIndex,G4int vacancyIndex);

  void LoadData( G4int Z);

  void PrintData();

private:

  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > idMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > probabilityMap;
  G4std::vector<G4int> nInitShells;
  G4int numberOfVacancies;
  
};

#endif


