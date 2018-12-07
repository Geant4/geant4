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
// File name:     G4EmDataHandler
//
// Author:        V. Ivanchenko 
// 
// Creation date: 16 August 2017
//
// Modifications: 
//
// -------------------------------------------------------------------
//
//    

#include "G4EmDataHandler.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4VEmProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::G4EmDataHandler(size_t n) : tLength(n)
{
  data.resize(n, nullptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::~G4EmDataHandler() 
{
  //std::cout << "G4EmDataHandler::~G4EmDataHandler " 
  //	    << tLength << "  "  << this << std::endl;
  for(size_t i=0; i<tLength; ++i) { CleanTable(i); }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t G4EmDataHandler::SetTable(G4PhysicsTable* ptr)
{
  data.push_back(ptr);
  ++tLength;
  return tLength-1; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsTable* G4EmDataHandler::MakeTable(size_t i)
{
  G4PhysicsTable* table = nullptr;
  // create new table only if index corresponds to the
  // position in the vector
  if(i <= tLength) {
    if(i < tLength) { table = data[i]; }
    table = G4PhysicsTableHelper::PreparePhysicsTable(table);
    if(i == tLength) {
      data.push_back(table);
      ++tLength;
    } else { data[i] = table; }
  }
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::CleanTable(size_t i)
{
  //std::cout << i << "  " << data[i] << std::endl;
  if(i < tLength && data[i]) {
    data[i]->clearAndDestroy();
    delete data[i];
    data[i] = nullptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::StorePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& fname, 
			   G4bool ascii)
{
  G4bool yes = true;
  if(data[idx]) {
    yes = data[idx]->StorePhysicsTable(fname, ascii);

    if ( yes ) {
      G4cout << "Physics table is stored for " 
	     << part->GetParticleName()
             << " <" << fname << "> " << G4endl;
    } else {
      G4cout << "Fail to store Physics Table for " 
             << part->GetParticleName()
             << " <" << fname << "> " << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::RetrievePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& fname, 
			   G4bool ascii)
{
  G4bool yes = 
    G4PhysicsTableHelper::RetrievePhysicsTable(data[idx], fname, ascii);
  G4EmParameters* param = G4EmParameters::Instance();
  if ( yes ) {
    if (0 < param->Verbose()) {
      G4cout << "Physics table " << idx << " for "  
	     << part->GetParticleName()
	     << " is retrieved from <" << fname << ">"
	     << G4endl;
    }
    if(param->Spline()) {
      G4PhysicsTable* table = data[idx];
      size_t n = table->length();
      for(size_t i=0; i<n; ++i) {
	if((*table)[i]) { (*table)[i]->SetSpline(true); } 
      }
    }
  } else if (1 < param->Verbose()) {
    G4cout << "Fail to retrieve physics table " << idx << " for " 
	   << part->GetParticleName() << " from <"
	   << fname << ">" << G4endl;
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::SetMasterProcess(const G4VEmProcess* ptr)
{
  masterProcess.push_back(ptr);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const G4VEmProcess* G4EmDataHandler::GetMasterProcess(size_t idx) const
{
  const G4VEmProcess* ptr = 
    (idx < masterProcess.size()) ? masterProcess[idx] : nullptr;
  return ptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
