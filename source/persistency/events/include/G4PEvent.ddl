// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PEvent.ddl,v 1.8 1999/11/28 21:54:17 morita Exp $
// GEANT4 tag $Name: geant4-01-01 $
//

#ifndef G4PEvent_h
#define G4PEvent_h 1

#include "globals.hh"
#include "G4PersistentTypes.hh"
#include "G4PersistentSchema.hh"

#include "G4PPrimaryVertex.hh"
#include "G4PHCofThisEvent.hh"
#include "G4PDCofThisEvent.hh"

class G4Event;

#include "HepODBMS/odbms/HepODBMS.h"


class G4PEvent 
  : public HepPersObj

{
  public:
      G4PEvent();
      G4PEvent(const G4Event *evt);
      G4PEvent(const G4Event *evt,
               HepRef(G4PHCofThisEvent) pHC,
               HepRef(G4PDCofThisEvent) pDC);
      void InitPEvent(const G4Event *evt,
                      HepRef(G4PHCofThisEvent) pHC,
                      HepRef(G4PDCofThisEvent) pDC);
      ~G4PEvent();

  private:
      // event ID
      G4Pint eventID;      

      // PrimaryVertex
      d_Ref<G4PPrimaryVertex> thePrimaryVertex;
      G4Pint numberOfPrimaryVertex;

      // HitsCollection
      d_Ref<G4PHCofThisEvent> HC;

      // DigiCollection
      d_Ref<G4PDCofThisEvent> DC;

      // TrajectoryContainer (not supported yet)
//      d_Ref<G4PTrajectoryContainer> trajectoryContainer;

  public:
      void SetEventID(const G4Event *evt);
      inline G4int GetEventID() const
      { return eventID; }

      inline void AddPrimaryVertex( HepRef(G4PPrimaryVertex) aPrimaryVertex)
      {
        if( thePrimaryVertex == NULL )
        { thePrimaryVertex = aPrimaryVertex; }
        else
        { thePrimaryVertex->SetNext( aPrimaryVertex ); }
        numberOfPrimaryVertex++;
      }

      inline G4int GetNumberOfPrimaryVertex() const
      { return numberOfPrimaryVertex; }

      inline HepRef(G4PPrimaryVertex) GetPrimaryVertex(G4int i=0)  const
      {
        if( i == 0 )
        { return thePrimaryVertex; }
        else if( i > 0 && i < numberOfPrimaryVertex )
        {
          HepRef(G4PPrimaryVertex) primaryVertex = thePrimaryVertex;
          for( int j=0; j<i; j++ )
          {
            if( primaryVertex == NULL ) return NULL;
            primaryVertex = primaryVertex->GetNext();
          }
          return primaryVertex;
        }
        else
        { return NULL; }
      }

      inline void SetHCofThisEvent( HepRef(G4PHCofThisEvent) value)
      { HC = value; }
      inline HepRef(G4PHCofThisEvent) GetHCofThisEvent()  const
      { return HC; }

      inline void SetDCofThisEvent( HepRef(G4PDCofThisEvent) value)
      { DC = value; }
      inline HepRef(G4PDCofThisEvent) GetDCofThisEvent()  const
      { return DC; }

      // inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      // { trajectoryContainer = value; }
      // inline G4TrajectoryContainer* GetTrajectoryContainer() const
      // { return trajectoryContainer; }

      G4Event* MakeTransientObject();

};

#endif

