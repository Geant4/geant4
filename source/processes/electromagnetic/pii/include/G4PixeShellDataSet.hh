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
// $Id: G4PixeShellDataSet.hh 70904 2013-06-07 10:34:25Z gcosmo $
//
// Author: Maria Grazia Pia (Maria.Grazia.Pia@cern.ch)
//
// History:
// -----------
// 31 Jul 2001   MGP        Created as G4EMShellDataSe
//  9 Mar 2008   MGP        Cleaned up unreadable code modified by former developer
//                          (Further clean-up needed) 
// 31 Jul 2008   MGP        Revised and renamed to G4PixeShellDataSet
//
// -------------------------------------------------------------------

// Class description:
// Shell data set 
// Applies a Composite design pattern for data library management

// -------------------------------------------------------------------

#ifndef  G4PIXESHELLDATASET_HH
#define  G4PIXESHELLDATASET_HH 1

#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>

#include "globals.hh"
#include "G4IDataSet.hh"

class G4IInterpolator;

class G4PixeShellDataSet : public G4IDataSet 
{ 
public:
  G4PixeShellDataSet(G4int Z, 
		     G4IInterpolator* algo, 
		     //const G4String& particleTpye="proton",
                     const G4String& modelK="ecpssr",
		     const G4String& modelL="ecpssr",
		     const G4String& modelM="ecpssr",
		     G4double eUnit=CLHEP::MeV, 
		     G4double dataUnit=CLHEP::barn);

  virtual ~G4PixeShellDataSet();
 
  virtual G4double FindValue(G4double energy, G4int componentId=0) const;
  
  virtual void PrintData(void) const;

  virtual const G4IDataSet* GetComponent(G4int componentId) const 
  { return components[componentId]; }

  virtual void AddComponent(G4IDataSet* dataSet) 
  { components.push_back(dataSet); }

  virtual size_t NumberOfComponents(void) const 
  { return components.size(); }

  virtual const G4DataVector& GetEnergies(G4int componentId) const 
  { return GetComponent(componentId)->GetEnergies(0); }

  virtual const G4DataVector& GetData(G4int componentId) const 
  { return GetComponent(componentId)->GetData(0); }

  virtual void SetEnergiesData(G4DataVector* energies, 
			       G4DataVector* data, 
			       G4int componentId);

  virtual G4bool LoadData(const G4String& fileName);
  virtual G4bool SaveData(const G4String& fileName) const;

  virtual G4double RandomSelect(G4int /*componentId = 0*/) const {return -1.;};
   
protected:

  G4double GetUnitEnergies() const { return unitEnergies; }
  G4double GetUnitData() const { return unitData; }
  const G4IInterpolator* GetAlgorithm() const { return algorithm; }
   
  void CleanUpComponents(void);

private:

  G4String FullFileName(const G4String& particleType,
			const G4String& subShell) const;
  G4int TranslateShell(const G4String& subShell) const;
  
  // Hide copy constructor and assignment operator 
  G4PixeShellDataSet();
  G4PixeShellDataSet(const G4PixeShellDataSet& copy);
  G4PixeShellDataSet& operator=(const G4PixeShellDataSet& right);

  std::vector<G4IDataSet*> components;          // Owned pointers

  G4int z;
  G4IInterpolator* algorithm;           // Owned pointer 
  // G4String particle;
  // G4String crossModelK;
  // G4String crossModelL;
  // G4String crossModelM;
  std::vector<G4String> crossModel;
  G4double unitEnergies;
  G4double unitData;
  std::vector<G4String> shellName;
  std::vector<G4String> subShellName;

};
#endif /* G4PIXESHELLDATASET_HH */
