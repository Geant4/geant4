// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PEvent.ddl,v 1.10 2000/12/05 14:40:02 morita Exp $
// GEANT4 tag $Name: geant4-03-01 $
//

// Class Description:
//
//      This is a class which represent a persistent event in Geant4.
//    A persistent event object is constructed by
//    G4PersistentEventMan::Store() or read back from a database
//    by G4PersistentEventMan::Retrieve().
//      G4PersistentEventMan::Store() passes a pointer of a transient
//    G4Event object.  The G4PEvent object construct itself by copying
//    the data member of G4Event, and by creating corresponting
//    persistent objects of primary vertex, particles, hits and digits.
//

#ifndef G4PEvent_h
#define G4PEvent_h 1

#include "G4Pglobals.hh"
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
  public: // with description
      G4PEvent();
      G4PEvent(const G4Event *evt);
      G4PEvent(const G4Event *evt,
               HepRef(G4PHCofThisEvent) pHC,
               HepRef(G4PDCofThisEvent) pDC);
      // Constructors.
      void InitPEvent(const G4Event *evt,
                      HepRef(G4PHCofThisEvent) pHC,
                      HepRef(G4PDCofThisEvent) pDC);
      // Actual constructor to be called by the above.
      ~G4PEvent();
      // Destructor.

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

  public: // with description
      void SetEventID(const G4Event *evt);
      inline G4int GetEventID() const
      { return eventID; }
      // Set and get event ID.

      inline void AddPrimaryVertex( HepRef(G4PPrimaryVertex) aPrimaryVertex)
      {
        if( thePrimaryVertex == NULL )
        { thePrimaryVertex = aPrimaryVertex; }
        else
        { thePrimaryVertex->SetNext( aPrimaryVertex ); }
        numberOfPrimaryVertex++;
      }
      // Adds another primary vertex to the vertex chain.

      inline G4int GetNumberOfPrimaryVertex() const
      { return numberOfPrimaryVertex; }
      // Get the number of primary vertecies in the chain.

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
      // returns a smart pointer of the i-th primary vertex.

      inline void SetHCofThisEvent( HepRef(G4PHCofThisEvent) value)
      { HC = value; }
      inline HepRef(G4PHCofThisEvent) GetHCofThisEvent()  const
      { return HC; }
      // Sets and gets a smart pointer of the hits collections of this event.

      inline void SetDCofThisEvent( HepRef(G4PDCofThisEvent) value)
      { DC = value; }
      inline HepRef(G4PDCofThisEvent) GetDCofThisEvent()  const
      { return DC; }
      // Sets and gets a smart pointer of the digits collections of this event.

      // inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      // { trajectoryContainer = value; }
      // inline G4TrajectoryContainer* GetTrajectoryContainer() const
      // { return trajectoryContainer; }

      G4Event* MakeTransientObject();
      // Construct a transient G4Event object from the data members of
      // G4PEvent.

};

#endif

