// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Event.hh,v 1.1 1999-01-07 16:06:32 gunter Exp $
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

      void Print() const;
      void Draw() const;

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
      { eventID =  i; };
      inline G4int GetEventID() const
      { return eventID; };
      inline void AddPrimaryVertex(G4PrimaryVertex* aPrimaryVertex)
      {
        if( thePrimaryVertex == NULL )
        { thePrimaryVertex = aPrimaryVertex; }
        else
        { thePrimaryVertex->SetNext( aPrimaryVertex ); }
        numberOfPrimaryVertex++;
      };
      inline G4int GetNumberOfPrimaryVertex() const
      { return numberOfPrimaryVertex; };
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
      };
      inline void SetHCofThisEvent(G4HCofThisEvent*value)
      { HC = value; };
      inline G4HCofThisEvent* GetHCofThisEvent()  const
      { return HC; };
      inline void SetDCofThisEvent(G4DCofThisEvent*value)
      { DC = value; };
      inline G4DCofThisEvent* GetDCofThisEvent()  const
      { return DC; };
      inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      { trajectoryContainer = value; };
      inline G4TrajectoryContainer* GetTrajectoryContainer() const
      { return trajectoryContainer; };
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

