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
// $Id: G4ParameterManager.cc 91116 2015-06-20 12:33:45Z ihrivnac $

// Author: Ivana Hrivnacova, 18/06/2013  (ivana@ipno.in2p3.fr)

#include "G4ParameterManager.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"

// mutex in a file scope

namespace {
  //Mutex to lock master manager when merging parameters
  G4Mutex mergeMutex = G4MUTEX_INITIALIZER;
}

G4ParameterManager* G4ParameterManager::fgMasterInstance = nullptr;
G4ThreadLocal G4ParameterManager* G4ParameterManager::fgInstance = nullptr;

//_____________________________________________________________________________
G4ParameterManager* G4ParameterManager::Instance()
{
  if ( fgInstance == nullptr ) {
    G4bool isMaster = ! G4Threading::IsWorkerThread();
    fgInstance = new G4ParameterManager(isMaster);
  }
  
  return fgInstance;
}    

//_____________________________________________________________________________
G4ParameterManager::G4ParameterManager(G4bool isMaster)
 : fState("Parameter", isMaster),
   fMap() 
{
  if ( ( isMaster && fgMasterInstance ) || ( fgInstance ) ) {
    G4ExceptionDescription description;
    description 
      << "      " 
      << "G4ParameterAnalysisManager already exists." 
      << "Cannot create another instance.";
    G4Exception("G4ParameterAnalysisManager::G4ParameterAnalysisManager()",
                "Analysis_F001", FatalException, description);
  }
  if ( isMaster ) fgMasterInstance = this;
  fgInstance = this;
}

//_____________________________________________________________________________
G4ParameterManager::~G4ParameterManager()
{
  // delete only parameters create by the mager itself
  for ( auto it : fParametersToDelete ) {
    delete it;
  }
}

//
// public methods
//

//_____________________________________________________________________________
void  G4ParameterManager::RegisterParameter(G4VParameter* parameter)
{
  fMap[parameter->GetName()] = parameter;
}

//_____________________________________________________________________________
G4VParameter*  
G4ParameterManager::GetParameter(const G4String& name, G4bool warn) const
{
  // get G4VParammeter from the map
  auto it = fMap.find(name);
  if ( it == fMap.end() ) {
    if ( warn) {
      G4ExceptionDescription description;
      description << "      " << "parameter " << name << " does not exist.";
      G4Exception("G4ParameterManager::GetParameter", 
                  "Analysis_W011", JustWarning, description);
    }
    return 0;
  }

  return it->second;
}

// //_____________________________________________________________________________
// G4bool G4ParameterManager::MergeImpl(tools::histo::hmpi* hmpi) 
// {
//   // tbd
//   return false;
// }

//_____________________________________________________________________________
void G4ParameterManager::Merge() 
{
  // Do nothing if  there are no parameters registered
  // or if master thread
  if ( (! fMap.size()) ||  (! G4Threading::IsWorkerThread()) ) return;

  // The manager on mastter must exist
  if ( ! fgMasterInstance ) {
    G4ExceptionDescription description;
    description 
      << "      " << "No master G4ParameterManager instance exists." 
      << G4endl 
      << "      " << "Parameters will not be merged.";
      G4Exception("G4ParameterManager::Merge()",
                "Analysis_W031", JustWarning, description);
    return;
  }

  // The worker manager just merges its parameters to the master
  // This operation needs a lock
  G4cout << "Go to merge parameters" << G4endl; 
  G4AutoLock lock(&mergeMutex);

  // the other manager has the map with identical keys
  auto it = fMap.begin();
  for ( auto itMaster : fgMasterInstance->fMap ) {
    // G4VParameter* parameter = (it++)->second;
    // G4VParameter* masterParameter = itMaster->second;
    // masterParameter->Merge(*(parameter));
    itMaster.second->Merge(*((it++)->second));
  }  
  lock.unlock();
}

//_____________________________________________________________________________
void G4ParameterManager::Reset()
{
// Reset histograms and profiles
  
  for ( auto it : fMap ) {
    it.second->Reset();
  }  
}  
 

