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
/// \file parallel/ThreadsafeScorers/src/TSRun.cc
/// \brief Implementation of the TSRun class
//
//
// $Id: TSRun.cc 93110 2015-11-05 08:37:42Z jmadsen $
//
//
/// TSRun contains three collections of hits maps: a thread-local hits map,
///     a global atomic hits map (implemented as a static since TSRun is
///     implemented as a thread-local instance), and a global "mutex" hits map
///     (also implemented as a static). The thread-local hits map is the
///     same as you will find in many other examples. The atomics hits map
///     is the purpose of this example. Code-wise, the implementation looks
///     extremely similar to the thread-local version with the 3 primary
///     exceptions: (1) construction - there should only be one instance so
///     it should be a static member variable or a pointer/reference to a
///     single instance elsewhere in the code (stored in ActionInitialization,
///     for instance); (2) It does not need to, nor should be, summed in
///     G4Run::Merge(); and (3) destruction -- it should only be cleared by
///     the master thread since there is only one instance.
/// A "mutex" hits map is also included as reference for checking the results
///     accumulated by the thread-local hits maps and atomic hits maps. The
///     differences w.r.t. this hits maps are computed in
///     TSRunAction::EndOfRunAction
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "TSRun.hh"

#include "G4SDManager.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4TAtomicHitsMap<G4double>*> TSRun::fAtomicRunMaps;

std::map<G4String, TSRun::MutexHitsMap_t> TSRun::fMutexRunMaps;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRun::TSRun(const G4String& mfd_name)
: G4Run()
{
    ConstructMFD(mfd_name);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TSRun::~TSRun()
{
    //--- Clear HitsMap for RUN
    for(unsigned i = 0; i < fRunMaps.size(); ++i)
        delete fRunMaps[i];

    if(!G4Threading::IsWorkerThread())
    {
        for(unsigned i = 0; i < fAtomicRunMaps.size(); ++i)
            delete fAtomicRunMaps[i];

        fAtomicRunMaps.clear();
        fMutexRunMaps.clear();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//    clear all data members.
void TSRun::ConstructMFD(const G4String& mfdName)
{

    G4SDManager* SDman = G4SDManager::GetSDMpointer();
    //=================================================
    //  Initalize RunMaps for accumulation.
    //  Get CollectionIDs for HitCollections.
    //=================================================
    G4MultiFunctionalDetector* mfd =
        (G4MultiFunctionalDetector*)(SDman->FindSensitiveDetector(mfdName));
    //
    if ( mfd )
    {
      //--- Loop over the registered primitive scorers.
      for (G4int icol = 0; icol < mfd->GetNumberOfPrimitives(); icol++){
        // Get Primitive Scorer object.
        G4VPrimitiveScorer* scorer = mfd->GetPrimitive(icol);
        // collection name and collectionID for HitsCollection,
        // where type of HitsCollection is G4THitsMap in case
        // of primitive scorer.
        // The collection name is given by <MFD name>/<Primitive
        // Scorer name>.
        G4String collectionName = scorer->GetName();
        G4String fullCollectionName = mfdName+"/"+collectionName;
        G4int    collectionID = SDman->GetCollectionID(fullCollectionName);
        //
        if ( collectionID >= 0 ){
          G4cout << "++ " << fullCollectionName<< " id " << collectionID
                 << G4endl;
          // Store obtained HitsCollection information into data members.
          // And, creates new G4THitsMap for accumulating quantities during RUN.
          fCollNames.push_back(fullCollectionName);
          fCollIDs.push_back(collectionID);
          fRunMaps.push_back(new G4THitsMap<G4double>(mfdName,
                                                      collectionName));
          if(!G4Threading::IsWorkerThread())
          {
            fAtomicRunMaps.push_back(new G4TAtomicHitsMap<G4double>
                                     (mfdName, collectionName));
            fMutexRunMaps[fCollNames[collectionID]].clear();
          }
        } else {
          G4cout << "** collection " << fullCollectionName << " not found. "
                 <<G4endl;
        }
      }
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a TSRun.
void TSRun::RecordEvent(const G4Event* aEvent)
{
    G4Run::RecordEvent(aEvent);

  //=============================
  // HitsCollection of This Event
  //============================
  G4HCofThisEvent* HCE = aEvent->GetHCofThisEvent();
  if (!HCE) return;

  for(unsigned i = 0; i < fCollIDs.size(); ++i)
  {
    G4int fCollID = fCollIDs.at(i);
    //=======================================================
    // Sum up HitsMap of this Event into HitsMap of this RUN
    //=======================================================
    G4THitsMap<G4double>* EvtMap = 0;
    if ( fCollID >= 0 )           // Collection is attached to HCE
      EvtMap = static_cast<G4THitsMap<G4double>*>(HCE->GetHC(fCollID));
    else
      G4cout <<" Error EvtMap Not Found " << G4endl;

    if ( EvtMap )
    {
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *fRunMaps[fCollID] += *EvtMap;
      // atomic run map
      *fAtomicRunMaps[fCollID] += *EvtMap;
      // mutex run map
      static G4Mutex mtx = G4MUTEX_INITIALIZER;
      {
        G4AutoLock lock(&mtx);
        for(const auto& itr : *EvtMap)
        {
          fMutexRunMaps[fCollNames[fCollID]][itr.first] += *itr.second;
        }
      }
      //----------------------------------------------------------------//
    }

  }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Merge hits map from threads
void TSRun::Merge(const G4Run* aTSRun)
{
    const TSRun* localTSRun = static_cast<const TSRun*>(aTSRun);

    for(unsigned i = 0; i < fRunMaps.size(); ++i)
      *fRunMaps[i] += *localTSRun->fRunMaps[i];

    G4Run::Merge(aTSRun);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Access HitsMap.
// by full description of collection name, that is
// <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4THitsMap<G4double>* TSRun::GetHitsMap(const G4String& collName) const
{
  for(unsigned i = 0; i < fCollNames.size(); ++i)
  {
    if(collName == fCollNames[i])
        return fRunMaps[i];
  }

  G4Exception("TSRun", collName.c_str(), JustWarning,
              "GetHitsMap failed to locate the requested HitsMap");
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Access AtomicsHitsMap.
// by full description of collection name, that is
// <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4TAtomicHitsMap<G4double>*
TSRun::GetAtomicHitsMap(const G4String& collName) const
{
  for(unsigned i = 0; i < fCollNames.size(); ++i)
  {
    if(collName == fCollNames[i])
        return fAtomicRunMaps[i];
  }

  G4Exception("TSRun", collName.c_str(), JustWarning,
              "GetHitsMap failed to locate the requested AtomicHitsMap");
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// Access AtomicsHitsMap.
// by full description of collection name, that is
// <MultiFunctional Detector Name>/<Primitive Scorer Name>
TSRun::MutexHitsMap_t*
TSRun::GetMutexHitsMap(const G4String& collName) const
{
  if(fMutexRunMaps.find(collName) != fMutexRunMaps.end())
      return &fMutexRunMaps[collName];

  G4Exception("TSRun", collName.c_str(), JustWarning,
              "GetHitsMap failed to locate the requested MutexHitsMap");
  return nullptr;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
