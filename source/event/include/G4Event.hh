// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Event.hh,v 1.2 1999-11-05 04:16:15 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#ifndef G4Event_h
#define G4Event_h 1

#include "globals.hh"
#include "G4Allocator.hh"
#include "G4PrimaryVertex.hh"
#include "G4HCofThisEvent.hh"
#include "G4DCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"

// class description:
//
//  This is the class which represents an event. This class is constructed and
// deleted by G4RunManager (or its derived class). When G4Event object is passed
// to G4EventManager, G4Event must have one or more primary vertexes and primary 
// particle(s) associated to the vertex(es) as an input of simulating an event.
// As the consequences of simulating an event, G4Event has trajectories, hits
// collections, and/or digi collections. 

class G4VHitsCollection;
class G4Event 
{
  public:
      G4Event();
      G4Event(G4int evID);
      ~G4Event();

      inline void *operator new(size_t);
      inline void operator delete(void* anEvent);

      int operator==(const G4Event &right) const;
      int operator!=(const G4Event &right) const;

  public: // with description
      void Print() const;
      // Print the event ID (starts with zero and increments by one) to G4cout.
      void Draw() const;
      // Invoke Draw() methods of all stored trajectories, hits, and digits.
      // For hits and digits, Draw() methods of the concrete classes must be
      // implemented. Otherwise nothing will be drawn.

  private:
      // event ID
      G4int eventID;      

      // PrimaryVertex
      G4PrimaryVertex* thePrimaryVertex;
      G4int numberOfPrimaryVertex;

      // HitsCollection
      G4HCofThisEvent* HC;

      // DigiCollection
      G4DCofThisEvent* DC;

      // TrajectoryContainer
      G4TrajectoryContainer * trajectoryContainer;

  public:
      inline void SetEventID(G4int i)
      { eventID =  i; }
      inline void SetHCofThisEvent(G4HCofThisEvent*value)
      { HC = value; }
      inline void SetDCofThisEvent(G4DCofThisEvent*value)
      { DC = value; }
      inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      { trajectoryContainer = value; }
  public: // with description
      inline G4int GetEventID() const
      { return eventID; }
      //  Returns the event ID
      inline void AddPrimaryVertex(G4PrimaryVertex* aPrimaryVertex)
      {
        if( thePrimaryVertex == NULL )
        { thePrimaryVertex = aPrimaryVertex; }
        else
        { thePrimaryVertex->SetNext( aPrimaryVertex ); }
        numberOfPrimaryVertex++;
      }
      //  This method sets a new primary vertex. This method must be invoked 
      // exclusively by G4VPrimaryGenerator concrete class.
      inline G4int GetNumberOfPrimaryVertex() const
      { return numberOfPrimaryVertex; }
      //  Returns number of primary vertexes the G4Event object has.
      inline G4PrimaryVertex* GetPrimaryVertex(G4int i=0)  const
      { 
        if( i == 0 )
        { return thePrimaryVertex; }
        else if( i > 0 && i < numberOfPrimaryVertex )
        {
          G4PrimaryVertex* primaryVertex = thePrimaryVertex;
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
      //  Returns i-th primary vertex of the event.
      inline G4HCofThisEvent* GetHCofThisEvent()  const
      { return HC; }
      inline G4DCofThisEvent* GetDCofThisEvent()  const
      { return DC; }
      inline G4TrajectoryContainer* GetTrajectoryContainer() const
      { return trajectoryContainer; }
      //  These three methods returns the pointers to the G4HCofThisEvent
      // (hits collections of this event), G4DCofThisEvent (digi collections
      // of this event), and G4TrajectoryContainer (trajectory coonainer),
      // respectively.
};

extern G4Allocator<G4Event> anEventAllocator;

inline void* G4Event::operator new(size_t)
{
  void* anEvent;
  anEvent = (void*)anEventAllocator.MallocSingle();
  return anEvent;
}

inline void G4Event::operator delete(void* anEvent)
{
  anEventAllocator.FreeSingle((G4Event*)anEvent);
}

#endif

