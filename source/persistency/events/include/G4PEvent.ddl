// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PEvent.ddl,v 2.2 1998/07/12 03:01:14 urbi Exp $
// GEANT4 tag $Name: geant4-00 $
//

#ifndef G4PEvent_h
#define G4PEvent_h 1

#include "globals.hh"
#include "G4PersistentTypes.hh"
#include "G4Event.hh"
#include "G4PPrimaryVertex.hh"

#include "HepODBMS/odbms/HepODBMS.h"


class G4PEvent 
  : public HepPersObj

{
  public:
      G4PEvent();
      G4PEvent(const G4Event *evt);
      ~G4PEvent();

  private:
      // event ID
      G4Pint eventID;      

      // PrimaryVertex
      ooRef(G4PPrimaryVertex) thePrimaryVertex;
      G4Pint numberOfPrimaryVertex;

      // HitsCollection
//      G4HCofThisEvent* HC;

      // TrajectoryContainer
//      G4TrajectoryContainer * trajectoryContainer;

  public:
      inline void SetEventID(const G4Event *evt)
      { eventID = evt->GetEventID(); };
      inline G4int GetEventID() const
      { return eventID; };

      // inline void AddPrimaryVertex()
      // {
      // primaryVertex should be added here
      // ;
      // };
      // inline G4int GetNumberOfPrimaryVertex() const
      // { return numberOfPrimaryVertex; };
      // inline G4PrimaryVertex* GetPrimaryVertex(G4int i=0)  const
      // { 
      // primaryVertex should be added here
      // ;
      // };
      // inline void SetHCofThisEvent()
      // { HC = value; };
      // inline G4HCofThisEvent* GetHCofThisEvent()  const
      // { return HC; };
      // inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      // { trajectoryContainer = value; };
      // inline G4TrajectoryContainer* GetTrajectoryContainer() const
      // { return trajectoryContainer; };
};

#endif

