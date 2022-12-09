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
// -------------------------------------------------------------------
//
// GEANT4 Class file
//
// File name:     G4HadDataHandler
//
// Author:        V. Ivanchenko 
// 
// Creation date: 05 July 2022
//
// Modifications: 
//
// -------------------------------------------------------------------
//
//    

#include "G4HadDataHandler.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadDataHandler::G4HadDataHandler(std::size_t n) : tLength(n)
{
  data.resize(n, nullptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4HadDataHandler::~G4HadDataHandler() 
{
  for(std::size_t i=0; i<tLength; ++i) {
    for(std::size_t j = i+1; j<tLength; ++j) {
      if(data[j] == data[i]) { data[j] = nullptr; }
    }
    CleanTable(i);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadDataHandler::AddTable(G4PhysicsTable* ptr)
{
  data.push_back(ptr);
  ++tLength;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadDataHandler::UpdateTable(G4PhysicsTable* ptr, std::size_t idx)
{
  // update table pointer but not delete previous
  if(idx < tLength) { 
    if(ptr != data[idx]) { data[idx] = ptr; }
  } else {
    G4cout << "### G4HadDataHandler::UpdateTable fail for idx=" << idx 
           << " length=" << tLength << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadDataHandler::CleanTable(std::size_t idx)
{
  if(idx < tLength && nullptr != data[idx]) {
    data[idx]->clearAndDestroy();
    delete data[idx];
    data[idx] = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4HadDataHandler::SetMasterProcess(const G4VProcess* ptr)
{
  masterProcess.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4VProcess* G4HadDataHandler::GetMasterProcess(const size_t idx) const
{
  return (idx < masterProcess.size()) ? masterProcess[idx] : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
