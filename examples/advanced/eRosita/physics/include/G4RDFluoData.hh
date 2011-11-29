//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
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

#ifndef G4RDFLUODATA_HH
#define G4RDFLUODATA_HH 1

#include "globals.hh"
#include <vector>
#include <map>

class G4RDFluoTransition;
class G4DataVector;

class G4RDFluoData
{
public:

  G4RDFluoData();

  ~G4RDFluoData();

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
  G4int StartShellId(G4int initIndex, G4int vacancyIndex) const;

  // Given the indexes of the starting and final shells for the 
  // transition, returns the transition energy
  G4double StartShellEnergy(G4int initIndex, G4int vacancyIndex) const;

  // Given the indexes of the starting and final shells for the 
  // transition, returns the probability of this transition
  G4double StartShellProb(G4int initIndex, G4int vacancyIndex) const;

  void LoadData(G4int Z);

  void PrintData();

  //void BuildFluoTransitionTable();

  //std::vector<G4RDFluoTransition*> GetFluoTransitions(G4int Z);
  //G4RDFluoTransition GetFluoTransition(G4int Z, G4int shellId);

private:

  // Hide copy constructor and assignment operator 
  G4RDFluoData& operator=(const G4RDFluoData& right);
  G4RDFluoData(const G4RDFluoData&);

  std::map<G4int,G4DataVector*,std::less<G4int> > idMap;
  std::map<G4int,G4DataVector*,std::less<G4int> > energyMap;
  std::map<G4int,G4DataVector*,std::less<G4int> > probabilityMap;
  std::vector<G4int> nInitShells;
  G4int numberOfVacancies;
  std::map<G4int,std::vector<G4RDFluoTransition*>,std::less<G4int> > fluoTransitionTable;  
};

#endif


