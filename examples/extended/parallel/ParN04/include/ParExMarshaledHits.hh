//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: ParExMarshaledHits.hh,v 1.1 2002-03-05 15:22:13 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// --------------------------------------------------------------------
//                   Parallel Library for Geant4
//
//             Gene Cooperman <gene@ccs.neu.edu>, 2001
// --------------------------------------------------------------------

#ifndef ParExMarshaledHits_h
#define ParExMarshaledHits_h 1

#include "G4THitsCollection.hh"
#include "ParMarshaledObj.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

// Calling sequence:
// MarshaledHCofThisEvent::MarshaledHCofThisEvent()
//    -> MarshaledHCofThisEvent::MarshalHitsCollection
//                                       (G4VHitsCollection * aBaseHC, G4String HCname)
//    -> TMarshaledHitsCollection<T>::Marshal(G4VHitsCollection * aBaseHC)
//    -> TMarshaledHitsCollection<T>::MarshalHit(T *aHit)
// MarshaledHCofThisEvent::UnmarshalSlaveHCofThisEvent()
//    -> TMarshaledHitsCollection<T>::UnmarshalSlaveHitsCollection
//					(G4String HCname, G4String SDname)
//    -> TMarshaledHitsCollection<T>::Unmarshal(G4THitsCollection<T> * aHC)
//    -> TMarshaledHitsCollection<T>::UnmarshalHit(T * aHit)

// class T derived from G4VHit
template <class T> class TMarshaledHitsCollection
  : public MarshaledObj
{
  public:
    inline void Marshal(G4VHitsCollection * aBaseHC); // Marshal hits of aBaseHC

    inline void UnmarshalSlaveHitsCollection(G4String & HCname, G4String & SDname);
    inline void Unmarshal(G4THitsCollection<T> * aHC);

    // By default, we a hit of type T as a void pointer
    inline void MarshalHit( T *aHit )
    { MarshaledObj::Marshal((void *)aHit, sizeof(*aHit)); }

    // Unmarshal hits into HC for T.  By default assume it was marshaled
    //   as a void pointer (ObjPtr)
    inline void UnmarshalHit(T * aHit)
    {
      if (isNewHitsCollection)
        *aHit = *(T *)MarshaledObj::UnmarshalAsSharedObjPtr();
      else
        G4Exception(
          "Event Level Parallelism:  Every hit collection sent to master"
           " should be new.\n  Current hit collection, "
           + collectionName + ", is not new.\n"
        );
    }

  protected:
    static G4bool isNewHitsCollection;
    static G4String collectionName;
};

template <class T>
G4bool TMarshaledHitsCollection<T>::isNewHitsCollection = true;

template <class T>
G4String TMarshaledHitsCollection<T>::collectionName = "<not initialized>";

template <class T>
inline void TMarshaledHitsCollection<T>::Marshal(G4VHitsCollection * aBaseHC)
{ G4THitsCollection<T> *aHC = (G4THitsCollection<T> *)aBaseHC;
  G4int numHits = aHC->entries();
  MarshaledObj::Marshal( numHits );
#ifdef TOPC_DEBUG
  G4cout << "Marshaled numHits: " << numHits << G4endl;
#endif
  for (int j = 0; j < numHits; j++) {
    T *aHit = (*aHC)[j];
    MarshalHit( aHit );
  }
}

template <class T>
inline void TMarshaledHitsCollection<T>::Unmarshal
				(G4THitsCollection<T> * aHC)
{
  G4int numHits;
  MarshaledObj::Unmarshal( numHits );
#ifdef TOPC_DEBUG
  G4cout << "numHits: " << numHits << G4endl;
#endif
  for (int j = 0; j < numHits; j++) {
    T *aHit = new T;
    UnmarshalHit( aHit );
    aHC->insert( aHit );
  }
}

template <class T>
inline void TMarshaledHitsCollection<T>::UnmarshalSlaveHitsCollection
				(G4String & HCname, G4String & SDname)
{
  G4HCofThisEvent *HCE = G4RunManager::GetRunManager()->
                           GetCurrentEvent()->
                           GetHCofThisEvent();
  // Find or create HitsCollection, aHC, with name HCname
  G4THitsCollection<T> *aHC;
  // HCname is unique.
  collectionName = HCname;
  G4int HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  if (HCID < 0)
  { aHC = new G4THitsCollection<T>(SDname, HCname);
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
    HCE->AddHitsCollection( HCID, aHC );
    isNewHitsCollection = true;
  }
  else {
    aHC = (G4THitsCollection<T> *)HCE->GetHC(HCID);
    isNewHitsCollection = ( aHC->entries() == 0 ? true : false );
  }

  Unmarshal(aHC);
};

#endif
