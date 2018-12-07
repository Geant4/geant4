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
// SUSANNA: This class is based on the RE02 extended example
//

#include "ML2Run.hh"
#include "G4SDManager.hh"

#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"

ML2Run::ML2Run(const std::vector<G4String> mfdName) : G4Run()
{
  G4SDManager* pSDman = G4SDManager::GetSDMpointer();
  //=================================================
  //  Initalize RunMaps for accumulation.
  //  Get CollectionIDs for HitCollections.
  //=================================================
  G4int nMfd = mfdName.size();
  for ( G4int idet = 0; idet < nMfd ; idet++){  // Loop for all MFD.
    G4String detName = mfdName[idet];
    //--- Seek and Obtain MFD objects from SDmanager.
    G4MultiFunctionalDetector* mfd =
      (G4MultiFunctionalDetector*)(pSDman->FindSensitiveDetector(detName));
    //
    if ( mfd ){
      //--- Loop over the registered primitive scorers.
      for (G4int icol = 0; icol < mfd->GetNumberOfPrimitives(); icol++){
        // Get Primitive Scorer object.
        G4VPrimitiveScorer* scorer=mfd->GetPrimitive(icol);
        // collection name and collectionID for HitsCollection,
        // where type of HitsCollection is G4THitsMap in case of primitive 
        // scorer.
        // The collection name is given by <MFD name>/<Primitive Scorer name>.
        G4String collectionName = scorer->GetName();
        G4String fullCollectionName = detName+"/"+collectionName;
        G4int    collectionID = pSDman->GetCollectionID(fullCollectionName);
        //
        if ( collectionID >= 0 ){
          G4cout << "++ "<<fullCollectionName<< " id " << collectionID
                 << G4endl;
          // Store obtained HitsCollection information into data members.
          // And, creates new G4THitsMap for accumulating quantities during RUN.
          fCollName.push_back(fullCollectionName);
          fCollID.push_back(collectionID);
          fRunMap.push_back(new G4THitsMap<G4double>(detName,collectionName));
        }else{
          G4cout << "** collection " << fullCollectionName << " not found. "
                 << G4endl;
        }
      }
    }
  }
}

//    clear all data members.
ML2Run::~ML2Run()
{
  //--- Clear HitsMap for RUN
  G4int nMap = fRunMap.size();
  for ( G4int i = 0; i < nMap; i++){
    if(fRunMap[i] ) fRunMap[i]->clear();
  }
  fCollName.clear();
  fCollID.clear();
  fRunMap.clear();
}

//
//  RecordEvent is called at end of event.
//  For scoring purpose, the resultant quantity in a event,
//  is accumulated during a Run.
void ML2Run::RecordEvent(const G4Event* aEvent)
{
  numberOfEvent++;  // This is an original line.

  //=============================
  // HitsCollection of This Event
  //============================
  G4HCofThisEvent* pHCE = aEvent->GetHCofThisEvent();
  if (!pHCE) return;

  //=======================================================
  // Sum up HitsMap of this Event  into HitsMap of this RUN
  //=======================================================
  G4int nCol = fCollID.size();
  for ( G4int i = 0; i < nCol ; i++ ){  // Loop over HitsCollection
    G4THitsMap<G4double>* evtMap=0;
    if ( fCollID[i] >= 0 ){           // Collection is attached to pHCE
      evtMap = (G4THitsMap<G4double>*)(pHCE->GetHC(fCollID[i]));
    }else{
      G4cout <<" Error evtMap Not Found "<< i << G4endl;
    }
    if ( evtMap )  {
      //=== Sum up HitsMap of this event to HitsMap of RUN.===
      *fRunMap[i] += *evtMap;
      //======================================================
    }
  }
}

void ML2Run::Merge(const G4Run * aRun) 
{
  const ML2Run * localRun = static_cast<const ML2Run *>(aRun);
  
  //=======================================================
  // Merge HitsMap of working threads
  //=======================================================
  G4int nCol = localRun->fCollID.size();
  for ( G4int i = 0; i < nCol ; i++ ){  // Loop over HitsCollection
    if ( localRun->fCollID[i] >= 0 ){
      *fRunMap[i] += *localRun->fRunMap[i];
    }
  }
  
  G4Run::Merge(aRun);
}

//  Access method for HitsMap of the RUN
//
//-----
// Access HitsMap.
//  By  MultiFunctionalDetector name and Collection Name.
G4THitsMap<G4double>* ML2Run::GetHitsMap(const G4String& detName,
                                         const G4String& colName){
    G4String fullName = detName+"/"+colName;
    return GetHitsMap(fullName);
}

// Access HitsMap.
//  By full description of collection name, that is
//    <MultiFunctional Detector Name>/<Primitive Scorer Name>
G4THitsMap<G4double>* ML2Run::GetHitsMap(const G4String& fullName){
    G4int nCol = fCollName.size();
    for ( G4int i = 0; i < nCol; i++){
        if ( fCollName[i] == fullName ){
            return fRunMap[i];
        }
    }
    return NULL;
}

void ML2Run::DumpAllScorer(){

  // - Number of HitsMap in this RUN.
  G4int n = GetNumberOfHitsMap();
  // - GetHitsMap and dump values.
  for ( G4int i = 0; i < n ; i++ ){
    G4THitsMap<G4double>* runMap =GetHitsMap(i);
    if ( runMap ) {
      G4cout << " PrimitiveScorer RUN " 
             << runMap->GetSDname() <<","<< runMap->GetName() << G4endl;
      G4cout << " Number of entries " << runMap->entries() << G4endl;
    }
  }
}


