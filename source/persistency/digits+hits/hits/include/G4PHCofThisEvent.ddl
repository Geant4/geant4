// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PHCofThisEvent.ddl,v 1.8 1999/12/02 16:10:21 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4PHCofThisEvent_h
#define G4PHCofThisEvent_h 1

#include "G4PersistentSchema.hh"

#include "G4PVHitsCollection.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PHCofThisEvent 
 : public HepPersObj
{
  public:
      G4PHCofThisEvent();
      ~G4PHCofThisEvent();

      void AddHitsCollection(G4int HCID, HepRef(G4PVHitsCollection) aHC);

  private:
      d_Varray< d_Ref<G4PVHitsCollection> > HC;

  public:
      inline HepRef(G4PVHitsCollection) GetHC(G4int i)
      { return HC[i]; }
      inline G4int GetCapacity()
      { return HC.size(); }
      inline G4int GetNumberOfCollections()
      {
        G4int n = 0;
        for(int i=0;i<HC.size();i++)
        {
          if( HC[i] != NULL ) n++;
        }
        return n;
      }
};

#endif

