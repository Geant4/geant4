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
// $Id: G4NistMaterialManager.hh,v 1.1 2005-02-11 17:30:24 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
// -------------------------------------------------------------------
//
// GEANT4 Class header file
//
// File name:     G4NistMaterialManager
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
// A utility static class, responsable for G4Element and G4Material access
//

// -------------------------------------------------------------------
//

#ifndef G4NistMaterialManager_h
#define G4NistMaterialManager_h 1

#include <vector>
#include "globals.hh"
#include "G4Material.hh"
#include "G4NistElementBuilder.hh"

///class G4NistMaterialMessenger;
class G4NistMaterialBuilder;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

class G4NistMaterialManager
{

public:

   G4NistMaterialManager(G4int verb=0);
  ~G4NistMaterialManager();

 // Elements
 //
  void RegisterElement  (const G4Element*);
  void DeRegisterElement(const G4Element*);

  const G4Element* GetElement(size_t index);
  
  // Find or build G4Element by atomic number
  const G4Element* FindOrBuildElement(G4int Z, G4bool isotopes=true);
  
 // Find or build G4Element by symbol
  const G4Element* FindOrBuildElement(const G4String& symb,
                                      G4bool isotopes=true);

  size_t   GetNumberOfElements() {return nElements;};
  G4int    GetZ (const G4String& symb);
  G4double GetIsotopeMass (G4int Z, G4int N);

  void PrintElement (const G4String&);
  void PrintElement (G4int Z);

  void PrintIsotopes (const G4String&);
  void PrintIsotopes (G4int Z);

  void ListElements (G4bool isotopes = true);

 // Materials
 //
  void RegisterMaterial  (G4Material*);
  void DeRegisterMaterial(G4Material*);

  G4Material* GetMaterial(size_t index);
  
  // Find or build a G4Material by name, from the dataBase
  //
  G4Material* FindOrBuildMaterial(const G4String& name, G4bool isotopes=true);

  G4Material* ConstructNewMaterial(const G4String& name,
                                      const std::vector<G4int>& Z,
                                      const std::vector<G4double>& atomFraction,
				      G4double dens, G4bool isotopes=true);

  G4Material* ConstructNewMaterial(const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& atomFraction,
				      G4double dens, G4bool isotopes=true);

  size_t GetNumberOfMaterials() {return nMaterials;};
  
  void SetVerbose(G4int);
  void PrintMaterial(const G4String&);
  void ListNistSimpleMaterials();
  void ListNistCompoundMaterials();
  void ListHepMaterials();

private:

  std::vector<const G4Element*>  elements;
  std::vector<G4Material*>       materials;
  
  size_t   nElements;
  size_t   nMaterials;
  
  G4int    verbose;
  
  G4NistElementBuilder*    elmBuilder;
  G4NistMaterialBuilder*   matBuilder;
  ///G4NistMaterialMessenger* messenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline const G4Element* G4NistMaterialManager::GetElement(size_t index)
{
  const G4Element* elm = 0;
  if(index < nElements) elm = elements[index];
  return elm;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
const G4Element* G4NistMaterialManager::FindOrBuildElement(G4int Z,
                                                           G4bool isotopes)
{
  return elmBuilder->FindOrBuildElement(Z, isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline 
const G4Element* G4NistMaterialManager::FindOrBuildElement(const G4String& symb,
                                                             G4bool isotopes)
{
  return elmBuilder->FindOrBuildElement(symb, isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int G4NistMaterialManager::GetZ(const G4String& symb)
{
  return elmBuilder->GetZ(symb);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4double G4NistMaterialManager::GetIsotopeMass(G4int Z, G4int N)
{
  return elmBuilder->GetIsotopeMass(Z, N);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4Material* G4NistMaterialManager::GetMaterial(size_t index)
{
  G4Material* mat = 0;
  if(index < nMaterials) mat = materials[index];
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

