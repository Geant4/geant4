// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PDCofThisEvent.ddl,v 1.4 2000/12/15 08:04:13 morita Exp $
// GEANT4 tag $Name: geant4-03-00 $
//

// Class Description:
//   This is a persistent version of container class which stores
// the associations to the digi collection in one event.
//   User should check the existence of the instance of this class
// in his/her sensitive detector class with GetCurrentPDCofThisEvent()
// method of a singleton class G4PersistentDigitMan.  If it does
// not exist, use should create a new instance and register the smart
// pointer with SetCurrentPDCofThisEvent() method of G4PersistentDigitMan.
//   As user creates a new collection of digit, its smart pointer should
// be registered to this class via AddDigitsCollection().
//   To obtain a smart pointer of the i-th collection, use GetDC() method.
//

#ifndef G4PDCofThisEvent_h
#define G4PDCofThisEvent_h 1

#include "G4PersistentSchema.hh"

#include "G4PVDigitsCollection.hh"

#include "HepODBMS/odbms/HepODBMS.h"

class G4PDCofThisEvent 
 : public HepPersObj
{
  public: // with description
      G4PDCofThisEvent();
      // Constructor.
      ~G4PDCofThisEvent();
      // Destructor.
      void AddDigitsCollection(G4int DCID, HepRef(G4PVDigitsCollection) aDC);
      // Store a smart pointer of the collection aDC as DCID-th entry.

  private:
      d_Varray< d_Ref<G4PVDigitsCollection> > DC;

  public: // with description
      inline HepRef(G4PVDigitsCollection) GetDC(G4int i)
      { return DC[i]; }
      //  Returns a smart pointer of a digit collection at the index i.
      // Null will be returned if the particular collection is not stored
      // in the current event.

  public:
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

