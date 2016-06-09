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
// $Id: G4NistManager.cc,v 1.15 2007/12/11 13:32:08 gcosmo Exp $
// GEANT4 tag $Name: geant4-09-01 $
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
// 27.02.06 V.Ivanchneko add ConstructNewGasMaterial
// 18.04.06 V.Ivanchneko add combined creation of materials (NIST + user)
// 11.05.06 V.Ivanchneko add warning flag to FindMaterial method
// 26.07.07 V.Ivanchneko modify destructor to provide complete destruction
//                       of all elements and materials
// 27-07-07, improve destructor (V.Ivanchenko) 
// 28.07.07 V.Ivanchneko make simple methods inline
// 28.07.07 V.Ivanchneko simplify Print methods
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
#include "G4NistMessenger.hh"
#include "G4Isotope.hh"

G4NistManager* G4NistManager::instance = 0;
G4double G4NistManager::POWERZ13[256] = {0};
G4double G4NistManager::LOGA[256] = {0};

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
  matBuilder = new G4NistMaterialBuilder(elmBuilder,verbose);
  
  messenger  = new G4NistMessenger(this);  
  for(G4int i=1; i<256; i++) {
    G4double x = G4double(i);
    POWERZ13[i] = std::pow(x,1.0/3.0);
    LOGA[i] = std::log(x);
  }
  POWERZ13[0] = 1.0;
  LOGA[0]     = 0.0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4NistManager::~G4NistManager()
{
  //  G4cout << "NistManager: start material destruction" << G4endl;
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nmat = theMaterialTable->size();
  size_t i;
  for(i=0; i<nmat; i++) {
    if((*theMaterialTable)[i]) delete (*theMaterialTable)[i];
  }
  //  G4cout << "NistManager: start element destruction" << G4endl;
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  size_t nelm = theElementTable->size();
  for(i=0; i<nelm; i++) {
    if((*theElementTable)[i]) delete (*theElementTable)[i];
  }
  //  G4cout << "NistManager: start isotope destruction" << G4endl;
  const G4IsotopeTable* theIsotopeTable = G4Isotope::GetIsotopeTable();
  size_t niso = theIsotopeTable->size();
  for(i=0; i<niso; i++) {
    if((*theIsotopeTable)[i]) delete (*theIsotopeTable)[i];
  }
  //  G4cout << "NistManager: end isotope destruction" << G4endl;
  delete messenger;
  delete matBuilder;
  delete elmBuilder;  
  // G4cout << "NistManager: end destruction" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintElement(const G4String& symbol)
{
  if (symbol == "all") elmBuilder->PrintElement(0);
  else                 elmBuilder->PrintElement(elmBuilder->GetZ(symbol));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintG4Element(const G4String& name)
{
  const G4ElementTable* theElementTable = G4Element::GetElementTable();
  size_t nelm = theElementTable->size();
  for(size_t i=0; i<nelm; i++) {
    G4Element* elm = (*theElementTable)[i];
    if ( name == elm->GetName() || "all" == name) {
      G4cout << *elm << G4endl;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4NistManager::PrintG4Material(const G4String& name)
{
  const G4MaterialTable* theMaterialTable = G4Material::GetMaterialTable();
  size_t nmat = theMaterialTable->size();
  for(size_t i=0; i<nmat; i++) {
    G4Material* mat = (*theMaterialTable)[i];
    if ( name == mat->GetName() || "all" == name) {
      G4cout << *mat << G4endl;
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

G4double G4NistManager::GetZ13(G4double Z)
{
  G4int iz = G4int(Z);
  G4double x = (Z - G4double(iz))/(3.0*Z);
  if(iz > 255) iz = 255;
  else if(iz < 0) iz = 0;
  return POWERZ13[iz]*(1.0 + x);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double G4NistManager::GetLOGA(G4double A)
{
  G4int ia = G4int(A);
  G4double x = (A - G4double(ia))/A;
  if(ia > 255) ia = 255;
  else if(ia < 0) ia = 0;
  return LOGA[ia] + x;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
