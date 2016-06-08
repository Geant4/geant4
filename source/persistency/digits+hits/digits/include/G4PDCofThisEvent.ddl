// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PDCofThisEvent.ddl,v 1.3 1999/12/02 16:10:19 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4PDCofThisEvent_h
#define G4PDCofThisEvent_h 1

#include "G4PersistentSchema.hh"

#include "G4PVDigitsCollection.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PDCofThisEvent 
 : public HepPersObj
{
  public:
      G4PDCofThisEvent();
      ~G4PDCofThisEvent();

      void AddDigitsCollection(G4int DCID, HepRef(G4PVDigitsCollection) aDC);

  private:
      d_Varray< d_Ref<G4PVDigitsCollection> > DC;

  public:
      inline HepRef(G4PVDigitsCollection) GetDC(G4int i)
      { return DC[i]; }
      inline G4int GetCapacity()
      { return DC.size(); }
      inline G4int GetNumberOfCollections()
      {
        G4int n = 0;
        for(int i=0;i<DC.size();i++)
        {
          if(! DC[i] ) n++;
        }
        return n;
      }
};

#endif

