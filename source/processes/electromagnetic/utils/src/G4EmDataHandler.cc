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
#include "G4EmDataRegistry.hh"
#include "G4ParticleDefinition.hh"
#include "G4EmParameters.hh"
#include "G4PhysicsTableHelper.hh"
#include "G4VEmProcess.hh"
#include "G4VEnergyLossProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::G4EmDataHandler(std::size_t n, const G4String& nam)
  : tLength(n), fName(nam)
{
  data.resize(n, nullptr);
  fMaxXS = new std::vector<G4double>;
  fXSpeaks = new std::vector<G4TwoPeaksXS*>;
  fRegistry = G4EmDataRegistry::Instance();
  fRegistry->Register(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4EmDataHandler::~G4EmDataHandler() 
{
  if (!fUseBaseParticleTable) {
    for (std::size_t i=0; i<tLength; ++i) {
      CleanTable(i);
    }
    delete fMaxXS;
    delete fXSpeaks;
  }
  if (!fElemSelectors.empty()) {
    for (auto const & ptr : fElemSelectors) {
      if (nullptr != ptr) {
	for (auto const & p : *ptr) { delete p; }
      }
    }
  }
  fRegistry->DeRegister(this);  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::size_t G4EmDataHandler::SetTable(G4PhysicsTable* ptr)
{
  for (std::size_t i=0; i<tLength; ++i) {
    if (ptr == data[i]) { return i; }
  }
  data.push_back(ptr);
  ++tLength;
  return tLength-1; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::UpdateTable(G4PhysicsTable* ptr, std::size_t idx)
{
  // update table pointer but not delete previous
  if (idx < tLength) { 
    if (ptr != data[idx]) {
      data[idx] = ptr;
    }
    data[idx] = G4PhysicsTableHelper::PreparePhysicsTable(data[idx]);
  } else {
    G4cout << "### G4EmDataHandler::UpdateTable fail for idx=" << idx 
           << " length=" << tLength << G4endl;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::SaveTable(G4PhysicsTable* ptr, std::size_t idx)
{
  if (idx < tLength) { 
    data[idx] = ptr;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsTable* G4EmDataHandler::MakeTable(std::size_t i)
{
  std::size_t idx = i; 
  if (idx >= tLength) { 
    data.push_back(nullptr);
    idx = tLength;
    ++tLength;
  }
  data[idx] = G4PhysicsTableHelper::PreparePhysicsTable(data[idx]);
  return data[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4PhysicsTable* G4EmDataHandler::MakeTable(G4PhysicsTable* ptr, std::size_t i)
{
  std::size_t idx = i; 
  // create new table only if index corresponds to the
  // position in the vector
  if (idx < tLength) { 
    if (ptr != data[idx]) { 
      CleanTable(idx);
      data[idx] = ptr;
    }
  } else {
    data.push_back(ptr);
    idx = tLength;
    ++tLength;
  }
  data[idx] = G4PhysicsTableHelper::PreparePhysicsTable(data[idx]);
  return data[idx];
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::CleanTable(std::size_t i)
{
  if (i < tLength && nullptr != data[i]) {
    auto ptr = data[i];
    ptr->clearAndDestroy();
    delete ptr;
    for (std::size_t j=0; j<tLength; ++j) {
      if (ptr == data[j]) { data[j] = nullptr; }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::StorePhysicsTable(std::size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& fname, 
			   G4bool ascii)
{
  G4bool yes = true;
  if(nullptr != data[idx]) {
    yes = data[idx]->StorePhysicsTable(fname, ascii);

    if ( yes ) {
      G4cout << "### Physics table is stored for " 
	     << part->GetParticleName()
             << " <" << fname << "> " << G4endl;
    } else {
      G4cout << "### Fail to store Physics Table for " 
             << part->GetParticleName()
             << " <" << fname << "> " << G4endl;
    }
  }
  return yes;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool G4EmDataHandler::RetrievePhysicsTable(std::size_t idx,
                           const G4ParticleDefinition* part,
			   const G4String& fname, 
			   G4bool ascii, G4bool spline)
{
  G4PhysicsTable* table = Table(idx);
  G4bool yes = G4PhysicsTableHelper::RetrievePhysicsTable(table, fname, ascii, spline);
  G4EmParameters* param = G4EmParameters::Instance();
  if ( yes ) {
    if (0 < param->Verbose()) {
      G4cout << "### Physics table " << idx << " for "  
	     << part->GetParticleName()
	     << " is retrieved from <" << fname << ">"
	     << G4endl;
    }
  } else if (1 < param->Verbose()) {
    G4cout << "### Fail to retrieve physics table " << idx << " for " 
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

const G4VEmProcess* G4EmDataHandler::GetMasterProcess(std::size_t idx) const
{
  return (idx < masterProcess.size()) ? masterProcess[idx] : nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void G4EmDataHandler::SetElementSelectors(std::vector<G4EmElementSelector*>* p,
					  std::size_t i)
{
  if (i < eLength) {
    if (fElemSelectors[i] != p) {
      delete fElemSelectors[i];
    }
    fElemSelectors[i] = p;
  } else {
    fElemSelectors.push_back(p);
    ++eLength;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
