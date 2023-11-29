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

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4AccumulableManager.hh"
#include "G4ThreadLocalSingleton.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging accumulables
  G4Mutex mergeMutex = G4MUTEX_INITIALIZER;
}

//_____________________________________________________________________________
G4AccumulableManager* G4AccumulableManager::Instance()
{
  static G4ThreadLocalSingleton<G4AccumulableManager> instance;
  return instance.Instance();
}

//_____________________________________________________________________________
G4AccumulableManager::G4AccumulableManager()
{
  if ( ! G4Threading::IsWorkerThread() ) fgMasterInstance = this;
}

//_____________________________________________________________________________
G4AccumulableManager::~G4AccumulableManager()
{
  // delete only accumulables create by the mager itself
  for ( auto it : fAccumulablesToDelete ) {
    delete it;
  }
}

//
// private methods
//

//_____________________________________________________________________________
G4String G4AccumulableManager::GenerateName() const
{
  G4String name = kBaseName;
  std::ostringstream os;
  os << fVector.size();
  name.append("_");
  name.append(os.str());
  return name;
}

//_____________________________________________________________________________
G4bool G4AccumulableManager::CheckName(const G4String& name, const G4String& where) const
{
  if ( fMap.find(name) == fMap.end() ) return true;

  G4ExceptionDescription description;
  description << "Name " << name << " is already used." << G4endl;
  description << "Parameter will be not created/registered.";
  G4String method("G4AccumulableManager::");
  method.append(where);
  G4Exception(method, "Analysis_W001", JustWarning, description);
  return false;
}

//
// public methods
//

//_____________________________________________________________________________
G4bool G4AccumulableManager::RegisterAccumulable(G4VAccumulable* accumulable)
{
  auto name = accumulable->GetName();

  // do not accept name if it is already used
  if ( ! CheckName(name, "RegisterAccumulable") ) return false;

  // generate name if empty
  if (name.length() == 0u) {
    name =  GenerateName();
    accumulable->fName = name;
  }

  fMap[name] = accumulable;
  fVector.push_back(accumulable);
  return true;
}

//_____________________________________________________________________________
G4VAccumulable*
G4AccumulableManager::GetAccumulable(const G4String& name, G4bool warn) const
{
  // get G4VParammeter from the map
  auto it = fMap.find(name);
  if ( it == fMap.end() ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "Accumulable " << name << " does not exist.";
      G4Exception("G4AccumulableManager::GetAccumulable",
                  "Analysis_W001", JustWarning, description);
    }
    return nullptr;
  }

  return it->second;
}

//_____________________________________________________________________________
G4VAccumulable*
G4AccumulableManager::GetAccumulable(G4int id, G4bool warn) const
{
  // get G4VParammeter from the vector
  if ( id < 0 || id >= G4int(fVector.size()) ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "Accumulable " << id << " does not exist.";
      G4Exception("G4AccumulableManager::GetAccumulable",
                  "Analysis_W001", JustWarning, description);
    }
    return nullptr;
  }

  return fVector[id];
}

//_____________________________________________________________________________
void G4AccumulableManager::Merge()
{
  // Do nothing if  there are no accumulables registered
  // or if master thread
  if ((fVector.size() == 0u) || (! G4Threading::IsWorkerThread())) return;

  // The manager on mastter must exist
  if (fgMasterInstance == nullptr) {
    G4ExceptionDescription description;
    description
      << "No master G4AccumulableManager instance exists." << G4endl
      << "Accumulables will not be merged.";
      G4Exception("G4AccumulableManager::Merge()",
                "Analysis_W001", JustWarning, description);
    return;
  }

  // The worker manager just merges its accumulables to the master
  // This operation needs a lock
  // G4cout << "Go to merge accumulables" << G4endl;
  G4AutoLock lock(&mergeMutex);

  // the other manager has the vector with the "same" accumulables
  auto it = fVector.begin();
  for ( auto itMaster : fgMasterInstance->fVector ) {
    // G4VAccumulable* masterAccumulable = itMaster;
    // G4VAccumulable* accumulable = *(it++);
    // masterAccumulable->Merge(*(accumulable));
    itMaster->Merge(*(*(it++)));
  }
  lock.unlock();
}

//_____________________________________________________________________________
void G4AccumulableManager::Reset()
{
// Reset histograms and profiles

  for ( auto it : fVector ) {
    it->Reset();
  }
}


