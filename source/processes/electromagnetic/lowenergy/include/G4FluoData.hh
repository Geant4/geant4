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
//  ?????
//
// -------------------------------------------------------------------

// Class description:
// Low Energy Electromagnetic Physics
// Fluorescence data set: shell identifiers, transition probabilities, 
// transition energies

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

  size_t NumberOfVacancies() const;

  G4int VacancyId(G4int vacancyIndex) const;

  size_t NumberOfTransitions(G4int vacancyIndex) const;

  const G4int StartShellId(G4int initIndex,G4int vacancyIndex);

  const G4double StartShellEnergy(G4int initIndex,G4int vacancyIndex);

  const G4double StartShellProb(G4int initIndex,G4int vacancyIndex);

  void LoadData(const G4int& Z);

  void PrintData();

private:

  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > idMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > energyMap;
  G4std::map<G4int,G4DataVector*,G4std::less<G4int> > probabilityMap;
  G4std::vector<G4int> nInitShells;
  G4int numberOfVacancies;
  
};

#endif


