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
//
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for cross sections per materials
//
// Author:      V.Ivanchenko 31.05.2018
//
// Modifications:
//
//----------------------------------------------------------------------------
//

#ifndef HadronXSDataTable_h
#define HadronXSDataTable_h 1

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "globals.hh"
#include "G4Element.hh"
#include "G4ElementVector.hh"
#include "G4PhysicsVector.hh"
#include "Randomize.hh"
#include <vector>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class G4Material;
class G4DynamicParticle;
class G4CrossSectionDataStore;

class G4HadElementSelector 
{
public:

  G4HadElementSelector(G4DynamicParticle*, G4CrossSectionDataStore*, 
		       const G4Material*, G4int bins, 
		       G4double emin, G4double emax, G4bool spline);

  ~G4HadElementSelector();

  void Dump();

  inline const G4Element* SelectRandomAtom(G4double e) const
  {
    const G4Element* element = (*theElementVector)[nElmMinusOne];
    if (nElmMinusOne > 0) {
      G4double x = G4UniformRand();
      for(G4int i=0; i<nElmMinusOne; ++i) {
	if (x <= xSections[i]->Value(e)) {
	  element = (*theElementVector)[i];
	  break;
	}
      }
    }
    return element;
  }

private:

  G4HadElementSelector(G4HadElementSelector &) = delete;
  G4HadElementSelector& operator=(const G4HadElementSelector &right) = delete;

  G4int nElmMinusOne;
  const G4ElementVector* theElementVector;
  std::vector<G4PhysicsVector*> xSections;
};

class G4HadronXSDataTable
{

public:

  explicit G4HadronXSDataTable();

  ~G4HadronXSDataTable();

  void Initialise(G4DynamicParticle*, G4CrossSectionDataStore*, 
		  G4int bins, G4double emin, G4double emax, 
		  G4bool spline);

  inline const G4PhysicsVector* HasData(size_t idx) const
  {
    return xsData[idx];
  };

  inline G4double GetCrossSection(G4double e, size_t idx) const
  {
    return xsData[idx]->Value(e);
  };

  inline const G4Element* SelectRandomAtom(G4double e, size_t idx) const
  {
    return elmSelectors[idx]->SelectRandomAtom(e);
  };

  void Dump();

private:

  // Assignment operator and copy constructor
  G4HadronXSDataTable & operator=
    (const G4HadronXSDataTable &right) = delete;
  G4HadronXSDataTable(const G4HadronXSDataTable&) = delete;
  
  std::vector<G4PhysicsVector*> xsData;
  std::vector<G4HadElementSelector*> elmSelectors;

  size_t nMaterials;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
 
