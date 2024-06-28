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
//---------------------------------------------------------------------------
//
// GEANT4 Class file
//
// Description: Data structure for registration of static cross sections components
//
// Author:      V.Ivanchenko 31.05.2018
//
// Modifications:
//
//----------------------------------------------------------------------------
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4HadronXSDataTable.hh"

G4HadronXSDataTable* G4HadronXSDataTable::sInstance = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronXSDataTable* G4HadronXSDataTable::Instance() {
  if ( sInstance == nullptr ) {
    static G4HadronXSDataTable theObject;
    sInstance = &theObject;
  }
  return sInstance;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronXSDataTable::G4HadronXSDataTable()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadronXSDataTable::~G4HadronXSDataTable()
{
  for (auto & ptr : fPiData) {
    delete ptr;
  }
  for (auto & ptr : fTable) {
    ptr->clearAndDestroy();
    delete ptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadronXSDataTable::AddPiData(std::vector<G4PiData*>* ptr)
{
  if (nullptr == ptr || ptr->empty()) { return; }
  for (auto & p : *ptr) {
    G4bool ok = true;
    for (auto & d : fPiData) {
      if (p == d) {
	ok = false;
	break;
      }
    }
    if (ok) { fPiData.push_back(p); }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadronXSDataTable::AddTable(G4PhysicsTable* ptr)
{
  if (nullptr != ptr) {
    for (auto & p : fTable) { if (p == ptr) { return; } }
    fTable.push_back(ptr);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
