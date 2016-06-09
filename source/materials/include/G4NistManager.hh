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
// $Id: G4NistManager.hh,v 1.3 2005/05/12 17:29:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4NistManager
//
// Author:        Vladimir Ivanchenko
//
// Creation date: 23.12.2004
//
// Modifications:
//
//
// Class Description:
//
// A utility static class
//

// -------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
// -------------------------------------------------------------------
//

#ifndef G4NistManager_h
#define G4NistManager_h 1

#include <vector>
#include "globals.hh"
#include "G4Material.hh"
#include "G4NistElementBuilder.hh"

class G4NistMaterialBuilder;
class G4NistMessenger;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4NistManager
{

public:

   static G4NistManager* Instance();
  ~G4NistManager();

 // Elements
 //
  void RegisterElement  (G4Element*);
  void DeRegisterElement(G4Element*);

  G4Element* GetElement(size_t index);
  
  // Find or build G4Element by atomic number
  G4Element* FindOrBuildElement(G4int Z, G4bool isotopes=true);
  
  // Find or build G4Element by symbol
  G4Element* FindOrBuildElement(const G4String& symb, G4bool isotopes=true);

  size_t   GetNumberOfElements() {return nElements;};
  G4int    GetZ (const G4String& symb);
  G4double GetIsotopeMass (G4int Z, G4int N);

  void PrintElement (const G4String&);
  void PrintElement (G4int Z);
    
  void PrintG4Element (const G4String&);  
  

 // Materials
 //
  void RegisterMaterial  (G4Material*);
  void DeRegisterMaterial(G4Material*);

  G4Material* GetMaterial(size_t index);
  
  // Find or build a G4Material by name, from the dataBase
  //
  G4Material* FindOrBuildMaterial(const G4String& name, G4bool isotopes=true);
  
  // construct a G4Material from scratch by atome count
  // 
  G4Material* ConstructNewMaterial(const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens, G4bool isotopes=true);
				      
  // construct a G4Material from scratch by fraction mass
  // 
  G4Material* ConstructNewMaterial(const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& weight,
				      G4double dens, G4bool isotopes=true);

  size_t GetNumberOfMaterials() {return nMaterials;};
  
  void SetVerbose(G4int);
  G4int GetVerbose();

  void ListMaterials(const G4String&);
  void PrintG4Material(const G4String&);

private:

  G4NistManager();
  static G4NistManager* instance;
  
  std::vector<G4Element*>   elements;
  std::vector<G4Material*>  materials;
  
  size_t   nElements;
  size_t   nMaterials;
  
  G4int    verbose;

  G4NistElementBuilder*    elmBuilder;
  G4NistMaterialBuilder*   matBuilder;
  G4NistMessenger*         messenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4Element* G4NistManager::GetElement(size_t index)
{
  G4Element* elm = 0;
  if(index < nElements) elm = elements[index];
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4Element* G4NistManager::FindOrBuildElement(G4int Z, G4bool isotopes)
{
  return elmBuilder->FindOrBuildElement(Z, isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4Element* G4NistManager::FindOrBuildElement(const G4String& symb,
                                                   G4bool isotopes)
{
  return elmBuilder->FindOrBuildElement(symb, isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4int G4NistManager::GetZ(const G4String& symb)
{
  return elmBuilder->GetZ(symb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistManager::GetIsotopeMass(G4int Z, G4int N)
{
  return elmBuilder->GetIsotopeMass(Z, N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4Material* G4NistManager::GetMaterial(size_t index)
{
  G4Material* mat = 0;
  if(index < nMaterials) mat = materials[index];
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline
G4int G4NistManager::GetVerbose()
{
  return verbose;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

