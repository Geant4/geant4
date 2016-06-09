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
// $Id: G4NistManager.cc,v 1.2 2005/05/12 17:29:08 vnivanch Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
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
// -------------------------------------------------------------------
//
// Class Description:
//
// Element data from the NIST DB on Atomic Weights and Isotope Compositions
// http://physics.nist.gov/PhysRefData/Compositions/index.html
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4NistManager.hh"

#include "G4NistMaterialBuilder.hh"
#include "G4NistMessenger.hh"

G4NistManager* G4NistManager::instance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4NistManager* G4NistManager::Instance()
{
  if (instance == 0) {
    static G4NistManager manager;
    instance = &manager;
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistManager::G4NistManager()
{
  nElements  = 0;
  nMaterials = 0;
  verbose    = 0;

  elmBuilder = new G4NistElementBuilder(verbose);
  matBuilder = new G4NistMaterialBuilder(this,elmBuilder,verbose);
  
  messenger  = new G4NistMessenger(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistManager::~G4NistManager()
{
  for (size_t i=0; i<nMaterials; i++) {
    if( materials[i] ) delete materials[i];
  }
  for (size_t j=0; j<nElements; j++) {
    if( elements[j] ) delete elements[j];
  }
  
  delete messenger;
  delete matBuilder;
  delete elmBuilder;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::RegisterElement(G4Element* elm)
{
  if(elm) {
    elements.push_back(elm);
    nElements++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::DeRegisterElement(G4Element* elm)
{
  for(size_t i=0; i<nElements; i++) {
    if(elm == elements[i]) {
      elements[i] = 0;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintElement(const G4String& symbol)
{
  if (symbol == "all") elmBuilder->PrintElement(0);
  else                 elmBuilder->PrintElement(elmBuilder->GetZ(symbol));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintElement(G4int Z)
{
  elmBuilder->PrintElement(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintG4Element(const G4String& name)
{
 for (size_t i=0; i<nElements; i++) {
  if ((name==(elements[i])->GetName()) || (name==(elements[i])->GetSymbol())) {
     G4cout << *(elements[i]) << G4endl;
     return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::RegisterMaterial(G4Material* mat)
{
  if(mat) {
    materials.push_back(mat);
    nMaterials++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::DeRegisterMaterial(G4Material* mat)
{
  for (size_t i=0; i<nMaterials; i++) {
    if (mat == materials[i]) { materials[i] = 0; return; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistManager::FindOrBuildMaterial(const G4String& name,
                                                       G4bool isotopes)
{
  if (verbose>1) G4cout << "G4NistManager::FindMaterial " << name 
                        << G4endl;
			
  G4Material* mat = 0;
  if (nMaterials > 0) {  
    for (size_t i=0; i<nMaterials; i++) {
       if (name == (materials[i])->GetName()) { mat = materials[i]; break; }
    }
  }
  if (!mat) mat = matBuilder->FindOrBuildMaterial(name, isotopes);  
  return mat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistManager::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4int>& nbAtoms,
				      G4double dens, G4bool isotopes)
{
  return matBuilder->ConstructNewMaterial(name,elm,nbAtoms,dens,isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistManager::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& w,
				      G4double dens, G4bool isotopes)
{
  return matBuilder->ConstructNewMaterial(name,elm,w,dens,isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::ListMaterials(const G4String& list)
{
  matBuilder->ListMaterials(list);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintG4Material(const G4String& name)
{
  for (size_t i=0; i<nMaterials; i++) {
    if (name == (materials[i])->GetName()) {
      G4cout << *(materials[i]) << G4endl;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::SetVerbose(G4int val)
{
  verbose = val;
  elmBuilder->SetVerbose(val);
  matBuilder->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
