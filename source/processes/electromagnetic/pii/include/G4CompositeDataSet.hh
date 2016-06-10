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
// $Id: G4CompositeDataSet.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created as G4CompositeEMDataSet
// 31 Jul 2008   MGP        Revised and renamed to G4CompositeDataSet
//
// -------------------------------------------------------------------

// Class description:
// Composite data set 
// Applies a Composite design pattern for data library management

// -------------------------------------------------------------------

#ifndef G4COMPOSITEDATASET_HH
#define G4COMPOSITEDATASET_HH 1

#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4IDataSet.hh"

class G4IInterpolator;

class G4CompositeDataSet : public G4IDataSet 
{
public:
  G4CompositeDataSet(G4IInterpolator* argAlgorithm, 
		     G4double eUnit=CLHEP::MeV, 
		     G4double dataUnit=CLHEP::barn, 
		     G4int zMin=1, 
		     G4int zMax=99); 

  virtual ~G4CompositeDataSet();
 
  virtual G4double FindValue(G4double x, G4int componentId=0) const;
  
  virtual void PrintData(void) const;

  virtual const G4IDataSet* GetComponent(G4int componentId) const 
  { return components[componentId]; }
  
  virtual void AddComponent(G4IDataSet* dataSet) 
  { components.push_back(dataSet); }
  
  virtual size_t NumberOfComponents() const 
  { return components.size(); }

  virtual const G4DataVector& GetEnergies(G4int componentId) const 
  { return GetComponent(componentId)->GetEnergies(0); }
  
  virtual const G4DataVector& GetData(G4int componentId) const 
  { return GetComponent(componentId)->GetData(0); }

  virtual void SetEnergiesData(G4DataVector* x, G4DataVector* data, G4int componentId);

  virtual G4bool LoadData(const G4String& fileName);
  virtual G4bool SaveData(const G4String& fileName) const;

  virtual G4double RandomSelect(G4int componentId) const;

private:

  void CleanUpComponents(void);
  
  // Hide copy constructor and assignment operator 
  G4CompositeDataSet();
  G4CompositeDataSet(const G4CompositeDataSet& copy);
  G4CompositeDataSet& operator=(const G4CompositeDataSet& right);

  std::vector<G4IDataSet*> components;    // Owned pointers

  G4IInterpolator* algorithm;           // Owned pointer 
  
  G4double unitEnergies;
  G4double unitData;

  G4int minZ;
  G4int maxZ;
};
#endif /* G4COMPOSITEDATASET_HH */









