// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4SolidStore.hh,v 1.2 1999-11-11 15:35:50 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// class G4SolidStore
//
// Container for all solids, with functionality derived from
// G4RWTPtrOrderedVector<T>. The class is `singleton', in that only
// one can exist, and access is facillitated via G4SolidStore::GetInstance()
//
// All solids should be registered with G4SolidStore, and removed on their
// destruction. Intended principally for UI browser. The underlying
// container initially has a capacity of 100.
//
// If much additional functionality is added, should consider containment
// instead of inheritance for G4RWTPtrOrderedVector<T>
//
// Class member functions:
//
// static void Register(G4G4VSolid* pSolid)
//   Add the solid to the collection
// static void DeRegister(G4G4VSolid* pSolid)
//   Remove the solid from the collection
// static G4SolidStore* GetInstance()
//   Get a ptr to the unique G4SolidStore, creting it if necessary
//
// Member functions:
//
// [as per RWTPtrOrderedvector]
//
// NOTE: Constructor is protected - creation and subsequent access is via
//       GetInstance
//
// Member data:
//
// static G4SolidStore*
//   Ptr to the single G4SolidStore
//
// History:
// 10.07.95 P.Kent Initial version

#ifndef G4VSOLIDSTORE_HH
#define G4VSOLIDSTORE_HH

#include "g4rw/tpordvec.h"

#include "G4VSolid.hh"

class G4SolidStore : public G4RWTPtrOrderedVector<G4VSolid>
{
  public:
    static void Register(G4VSolid* pSolid);
    static void DeRegister(G4VSolid* pSolid);
    static G4SolidStore* GetInstance();
    virtual ~G4SolidStore();
  protected:
    G4SolidStore();
  private:
    static G4SolidStore* fgInstance;
};

#endif
