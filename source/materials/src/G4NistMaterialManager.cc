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
// $Id: G4NistMaterialManager.cc,v 1.1 2005-02-11 17:30:25 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
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
// -------------------------------------------------------------------
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4NistMaterialManager.hh"

#include "G4NistMaterialBuilder.hh"
///#include "G4NistMaterialMessenger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistMaterialManager::G4NistMaterialManager(G4int vb)
{
  nElements  = 0;
  nMaterials = 0;
  verbose    = vb;

  elmBuilder = new G4NistElementBuilder(verbose);
  matBuilder = new G4NistMaterialBuilder(this,elmBuilder,verbose);
  
  ///messenger  = new G4NistMaterialMessenger();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistMaterialManager::~G4NistMaterialManager()
{
  for (size_t i=0; i<nMaterials; i++) {
    if( materials[i] ) delete materials[i];
  }
  for (size_t j=0; j<nElements; j++) {
    if( elements[j] ) delete elements[j];
  }
  ///delete messenger;
  delete matBuilder;
  delete elmBuilder;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::RegisterElement(const G4Element* elm)
{
  if(elm) {
    elements.push_back(elm);
    nElements++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::DeRegisterElement(const G4Element* elm)
{
  for(size_t i=0; i<nElements; i++) {
    if(elm == elements[i]) {
      elements[i] = 0;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::PrintElement(const G4String& name)
{
  G4int Z = elmBuilder->GetZ(name);
  elmBuilder->PrintElement(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::PrintElement(G4int Z)
{
  elmBuilder->PrintElement(Z);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::PrintIsotopes(const G4String& name)
{
  G4int Z = elmBuilder->GetZ(name);
  elmBuilder->PrintElement(Z, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::PrintIsotopes(G4int Z)
{
  elmBuilder->PrintElement(Z, true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::ListElements(G4bool isotopes)
{
  G4int nz = elmBuilder->GetMaxNumElements();
  for(G4int i=1; i<=nz; i++) {elmBuilder->PrintElement(i, isotopes);}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::RegisterMaterial(G4Material* mat)
{
  if(mat) {
    materials.push_back(mat);
    nMaterials++;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::DeRegisterMaterial(G4Material* mat)
{
  for (size_t i=0; i<nMaterials; i++) {
    if (mat == materials[i]) { materials[i] = 0; return; }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialManager::FindOrBuildMaterial(const G4String& name,
                                                       G4bool isotopes)
{
  if (verbose>1) G4cout << "G4NistMaterialManager::FindMaterial " << name 
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

G4Material* G4NistMaterialManager::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4int>& Z,
                                      const std::vector<G4double>& atomFraction,
				      G4double dens, G4bool isotopes)
{
  return matBuilder->ConstructNewMaterial(name,Z,atomFraction,dens,isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* G4NistMaterialManager::ConstructNewMaterial(
                                      const G4String& name,
                                      const std::vector<G4String>& elm,
                                      const std::vector<G4double>& atomFraction,
				      G4double dens, G4bool isotopes)
{
  return matBuilder->ConstructNewMaterial(name,elm,atomFraction,dens,isotopes);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::PrintMaterial(const G4String& name)
{
  for(size_t i=0; i<nMaterials; i++) {
    if(name == (materials[i])->GetName()) {
      G4cout << *(materials[i]) << G4endl;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::ListNistSimpleMaterials()
{
  matBuilder->ListNistSimpleMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::ListNistCompoundMaterials()
{
  matBuilder->ListNistCompoundMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::ListHepMaterials()
{
  matBuilder->ListHepMaterials();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistMaterialManager::SetVerbose(G4int val)
{
  verbose = val;
  elmBuilder->SetVerbose(val);
  matBuilder->SetVerbose(val);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
