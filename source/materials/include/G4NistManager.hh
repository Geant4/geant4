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
// $Id: G4NistManager.hh,v 1.9 2006/10/17 15:15:46 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-02 $
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
// 27.02.06 V.Ivanchneko add GetAtomicMassAmu and ConstructNewGasMaterial
// 11.05.06 V.Ivanchneko add warning flag to FindMaterial method
// 17.10.06 V.Ivanchneko add methods: GetAtomicMass, GetNistElementNames,
//                       GetNistMaterialNames
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

  size_t   GetNumberOfElements() const;
  G4int    GetZ(const G4String& symb) const;

  // Mass in amu
  G4double GetAtomicMassAmu(G4int Z) const;

  // Mass in Geant4 units
  G4double GetIsotopeMass(G4int Z, G4int N) const;
  G4double GetAtomicMass(G4int Z, G4int N) const;

  G4double GetIsotopeAbundance(G4int Z, G4int N) const;

  void PrintElement(const G4String&);
  void PrintElement(G4int Z);
    
  void PrintG4Element(const G4String&);  
  
  const std::vector<G4String>& GetNistElementNames() const;

  // Materials
  //
  void RegisterMaterial  (G4Material*);
  void DeRegisterMaterial(G4Material*);

  G4Material* GetMaterial(size_t index);
  
  // Find or build a G4Material by name, from the dataBase
  //
  G4Material* FindOrBuildMaterial(const G4String& name, G4bool isotopes=true,
				  G4bool warning=false);
  
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

  // construct a gas G4Material from scratch by atome count
  // 
  G4Material* ConstructNewGasMaterial(const G4String& name,
                                      const G4String& nameNist,
				      G4double temp, G4double pres, 
				      G4bool isotopes=true);

  size_t GetNumberOfMaterials() {return nMaterials;};
  
  void SetVerbose(G4int);
  G4int GetVerbose();

  void ListMaterials(const G4String&);
  void PrintG4Material(const G4String&);

  const std::vector<G4String>& GetNistMaterialNames() const;

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
size_t G4NistManager::GetNumberOfElements() const
{ 
  return nElements;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4int G4NistManager::GetZ(const G4String& symb) const
{
  return elmBuilder->GetZ(symb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistManager::GetAtomicMassAmu(G4int Z) const
{
  return elmBuilder->GetA(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistManager::GetIsotopeMass(G4int Z, G4int N) const
{
  return elmBuilder->GetIsotopeMass(Z, N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistManager::GetAtomicMass(G4int Z, G4int N) const
{
  return elmBuilder->GetAtomicMass(Z, N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
G4double G4NistManager::GetIsotopeAbundance(G4int Z, G4int N) const
{
  return elmBuilder->GetIsotopeAbundance(Z, N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
const std::vector<G4String>& G4NistManager::GetNistElementNames() const
{
  return elmBuilder->GetElementNames();
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

