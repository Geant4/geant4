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
// Author:  V.Ivanchenko 18.04.2024
//
//----------------------------------------------------------------------------

#include "G4EmDataRegistry.hh"
#include "G4AutoLock.hh"

G4EmDataRegistry* G4EmDataRegistry::instance = nullptr;

namespace
{
  G4Mutex theEmData = G4MUTEX_INITIALIZER;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....

G4EmDataRegistry* G4EmDataRegistry::Instance()
{
  if (instance == nullptr) {
    static G4EmDataRegistry manager;
    instance = &manager;
  }
  return instance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataRegistry::G4EmDataRegistry()
{
  fDataHandlers.reserve(50);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataRegistry::~G4EmDataRegistry()
{
  for (auto const & p : fDataHandlers) {
    delete p;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler*
G4EmDataRegistry::GetHandlerByName(const G4String& nam, std::size_t n)
{
  // handler already exist
  G4EmDataHandler* ptr = EmDataHandler(nam);
  if (nullptr != ptr) { return ptr; }

  // create a new handler
  G4AutoLock l(&theEmData);
  ptr = EmDataHandler(nam);
  if (nullptr == ptr) {
    ptr = new G4EmDataHandler(n, nam);
  }
  l.unlock();
  
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataRegistry::Register(G4EmDataHandler* ptr)
{
  if (nullptr == ptr) { return; }
  for (auto const & p : fDataHandlers) {
    if (p == ptr) { return; }
  }
  fDataHandlers.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......  

void G4EmDataRegistry::DeRegister(G4EmDataHandler* ptr)
{
  if (nullptr == ptr) { return; }
  std::size_t n = fDataHandlers.size();
  for (std::size_t i = 0; i < n; ++i) {
    if (fDataHandlers[i] == ptr) {
      fDataHandlers[i] = nullptr;
      return;
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler* G4EmDataRegistry::EmDataHandler(const G4String& nam)
{
  G4EmDataHandler* ptr = nullptr;
  if (fDataHandlers.empty()) { return ptr; }
  for (auto const & p : fDataHandlers) {
    if (p->GetName() == nam) {
      ptr = p;
      break;
    }
  }
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
