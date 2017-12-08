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
// $Id: G4EmDataHandler.cc 73844 2013-09-13 14:16:30Z vnivanch $
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

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::G4EmDataHandler(size_t n) : tLength(0)
{
  data.reserve(n);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::~G4EmDataHandler() 
{
  for(auto & table : data) { 
    if(table) {
      table->clearAndDestroy();
      delete table;
    }
  }
  data.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

size_t G4EmDataHandler::SetTable(G4PhysicsTable* ptr)
{
  data.push_back(ptr);
  ++tLength;
  return tLength-1; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsTable* G4EmDataHandler::MakeTable()
{
  G4PhysicsTable* table = nullptr;
  table = G4PhysicsTableHelper::PreparePhysicsTable(table);
  data.push_back(table);
  ++tLength;
  return table;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::StorePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& namet, 
			   const G4String& procname, 
			   G4bool ascii)
{
  G4bool yes = true;
  if(data[idx]) {
    yes = data[idx]->StorePhysicsTable(namet, ascii);

    if ( yes ) {
      G4cout << "Physics table is stored for " 
	     << part->GetParticleName()
             << " and process " << procname
             << " <" << namet << "> " << G4endl;
    } else {
      G4cout << "Fail to store Physics Table for " 
             << part->GetParticleName()
             << " and process " << procname
             << " <" << namet << "> " << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::RetrievePhysicsTable(size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& namet, 
			   const G4String& procname, 
			   G4bool ascii)
{
  G4bool yes = 
    G4PhysicsTableHelper::RetrievePhysicsTable(data[idx],
					       namet, ascii);
  G4EmParameters* param = G4EmParameters::Instance();
  if ( yes ) {
    if (0 < param->Verbose()) {
      G4cout << "Physics table for "  
	     << part->GetParticleName()
	     << " and " << procname
	     << " is retrieved from <" << namet << ">"
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
    G4cout << "Fail to retrieve physics table for " 
	   << part->GetParticleName()
	   << " and " << procname << " from <"
	   << namet << ">" << G4endl;
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
