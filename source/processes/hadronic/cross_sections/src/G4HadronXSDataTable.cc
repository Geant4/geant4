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
  for (std::size_t i = 0; i < fPiData.size(); ++i) {
    auto ptr = fPiData[i];
    for (std::size_t j = 0; j < ptr->size(); ++j) {
      auto p = (*ptr)[j];
      for (std::size_t k = i + 1; k < fPiData.size(); ++k) {
	auto qtr = fPiData[k];
	for (std::size_t l = 0; l < qtr->size(); ++l) {
	  if ((*qtr)[l] == p) { (*qtr)[l] = nullptr; }
	}
      }
      delete p;
      (*ptr)[j] = nullptr;
    }
    delete ptr;
  }
  fPiData.clear();
  for (auto const & ptr : fTable) {
    ptr->clearAndDestroy();
    delete ptr;
  }
  fTable.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadronXSDataTable::AddPiData(std::vector<G4PiData*>* ptr)
{
  if (nullptr == ptr || ptr->empty()) { return; }
  for (auto & d : fPiData) {
    if (ptr == d) { return; }
  }
  fPiData.push_back(ptr);
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
