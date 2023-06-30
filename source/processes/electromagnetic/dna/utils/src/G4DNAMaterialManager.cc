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
// 1/2/2023: Hoang: this file is used to check available DNA materials,
// and keeps DNA cross sections.

#include "G4DNAMaterialManager.hh"
#include "G4StateManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4NistManager.hh"

G4DNAMaterialManager *G4DNAMaterialManager::theInstance = nullptr;

namespace {
  G4Mutex MaterialMutex = G4MUTEX_INITIALIZER;
}

G4DNAMaterialManager *G4DNAMaterialManager::Instance() {
  if (nullptr == theInstance) {
    G4AutoLock l(&MaterialMutex);
    if (nullptr == theInstance) {
      static G4DNAMaterialManager manager;
      theInstance = &manager;
    }
    l.unlock();
  }
  return theInstance;
}

G4DNAMaterialManager::G4DNAMaterialManager() {
  G4NistManager::Instance();
  fStateManager = G4StateManager::GetStateManager();
}

G4VEmModel* G4DNAMaterialManager::GetModel(const DNAModelType& t){
  return fData[t];
}

void G4DNAMaterialManager::SetMasterDataModel(const DNAModelType& t, G4VEmModel* m){
  fData[t] = m;
}

G4bool G4DNAMaterialManager::IsLocked() const
{
  return (!G4Threading::IsMasterThread());
}
