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

//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Author:  V.Ivanchenko 26.08.2023
//
//----------------------------------------------------------------------------

#include "G4ElementDataRegistry.hh"

G4ElementDataRegistry* G4ElementDataRegistry::instance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4ElementDataRegistry* G4ElementDataRegistry::Instance()
{
  if (instance == nullptr) {
    static G4ElementDataRegistry manager;
    instance = &manager;
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElementDataRegistry::~G4ElementDataRegistry()
{
  for (auto const & ptr : elmdata) {
    delete ptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElementDataRegistry::G4ElementDataRegistry()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ElementDataRegistry::RegisterMe(G4ElementData* p)
{
  for (auto & ptr : elmdata) { if (ptr == p) { return; } } 

  elmdata.push_back(p);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4ElementDataRegistry::RemoveMe(G4ElementData* p)
{
  if (nullptr == p) { return; }
  for (std::size_t i=0; i<elmdata.size(); ++i) {
    if (p == elmdata[i]) {
      elmdata[i] = nullptr;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ElementData* G4ElementDataRegistry::GetElementDataByName(const G4String& nam)
{
  G4ElementData* ptr = nullptr;
  for (auto const & p : elmdata) {
    if (p->GetName() == nam) {
      ptr = p;
      break;
    }
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
