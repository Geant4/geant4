//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: G4Event.hh 94950 2016-01-07 11:53:14Z gcosmo $
//
//
// class description:
//
// This is the class which represents an event. A G4Event is constructed and
// deleted by G4RunManager (or its derived class). When a G4Event object is
// passed to G4EventManager, G4Event must have one or more primary verteces
// and primary particle(s) associated to the vertex(es) as an input of
// simulating an event.
// G4Event has trajectories, hits collections, and/or digi collections. 

#ifndef G4Event_h
#define G4Event_h 1

#include "globals.hh"
#include "evtdefs.hh"
#include "G4Allocator.hh"
#include "G4PrimaryVertex.hh"
#include "G4HCofThisEvent.hh"
#include "G4DCofThisEvent.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VUserEventInformation.hh"

class G4VHitsCollection;
class G4Event 
{
  public:
      G4Event();
      G4Event(G4int evID);
      ~G4Event();

      inline void *operator new(size_t);
      inline void operator delete(void* anEvent);

      G4int operator==(const G4Event &right) const;
      G4int operator!=(const G4Event &right) const;

  public: // with description
      void Print() const;
      // Print the event ID (starts with zero and increments by one) to G4cout.
      void Draw() const;
      // Invoke Draw() methods of all stored trajectories, hits, and digits.
      // For hits and digits, Draw() methods of the concrete classes must be
      // implemented. Otherwise nothing will be drawn.

  private:
      // Copy constructor and = operator must not be used.
      G4Event(const G4Event &) {;}
      G4Event& operator=(const G4Event &) { return *this; }

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

      // Boolean flag which shall be set to true if the event is aborted and 
      // thus the containing information is not to be used.
      G4bool eventAborted;

      // UserEventInformation (optional)
      G4VUserEventInformation* userInfo;

      // Initial random number engine status before primary particle generation
      G4String* randomNumberStatus;
      G4bool validRandomNumberStatus;

      // Initial random number engine status before event processing
      G4String* randomNumberStatusForProcessing;
      G4bool validRandomNumberStatusForProcessing;

      // Flag to keep the event until the end of run
      G4bool keepTheEvent;
      mutable G4int grips;

  public:
      inline void SetEventID(G4int i)
      { eventID =  i; }
      inline void SetHCofThisEvent(G4HCofThisEvent*value)
      { HC = value; }
      inline void SetDCofThisEvent(G4DCofThisEvent*value)
      { DC = value; }
      inline void SetTrajectoryContainer(G4TrajectoryContainer*value)
      { trajectoryContainer = value; }
      inline void SetEventAborted()
      { eventAborted = true; }
      inline void SetRandomNumberStatus(G4String& st)
      {
        randomNumberStatus = new G4String(st);
        validRandomNumberStatus = true;
      }
      inline void SetRandomNumberStatusForProcessing(G4String& st)
      {
        randomNumberStatusForProcessing = new G4String(st);
        validRandomNumberStatusForProcessing = true;
      }
      inline void KeepTheEvent(G4bool vl=true)
      { keepTheEvent = vl; }
      inline G4bool ToBeKept() const
      { return keepTheEvent; }
      inline void KeepForPostProcessing() const
      { grips++; }
      inline void PostProcessingFinished() const
      { grips--;
        if(grips<0)
        { G4Exception("G4Event::Release()","EVENT91001",FatalException,
                      "Number of grips becames negative. This cannot be correct."); }
      }
      inline G4int GetNumberOfGrips() const
      { return grips; }

  public: // with description
      inline G4int GetEventID() const
      { return eventID; }
      //  Returns the event ID
      inline void AddPrimaryVertex(G4PrimaryVertex* aPrimaryVertex)
      {
        if( thePrimaryVertex == nullptr )
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
          for( G4int j=0; j<i; j++ )
          {
            if( !primaryVertex ) return nullptr; 
            primaryVertex = primaryVertex->GetNext();
          }
          return primaryVertex;
        }
        else
        { return nullptr; }
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
      inline G4bool IsAborted() const { return eventAborted; }
      //  Return a boolean which indicates the event has been aborted and thus
      // it should not be used for analysis.
      inline void SetUserInformation(G4VUserEventInformation* anInfo) { userInfo = anInfo; }
      inline G4VUserEventInformation* GetUserInformation() const { return userInfo; }
      //  Set and Get method of G4VUserEventInformation
      inline const G4String& GetRandomNumberStatus() const 
      {
        if(!validRandomNumberStatus)
        { G4Exception(
              "G4Event::GetRandomNumberStatus","Event0701",JustWarning,
              "Random number status is not available for this event."); }
        return *randomNumberStatus;
      }
      inline const G4String& GetRandomNumberStatusForProcessing() const 
      {
        if(!validRandomNumberStatusForProcessing)
        { G4Exception(
              "G4Event::GetRandomNumberStatusForProcessing","Event0702",
              JustWarning,
              "Random number status is not available for this event."); }
        return *randomNumberStatusForProcessing;
      }
};

extern G4EVENT_DLL G4ThreadLocal G4Allocator<G4Event> *anEventAllocator;

inline void* G4Event::operator new(size_t)
{ 
  if (!anEventAllocator) anEventAllocator = new G4Allocator<G4Event>;
  return (void*)anEventAllocator->MallocSingle();
}

inline void G4Event::operator delete(void* anEvent)
{
  anEventAllocator->FreeSingle((G4Event*)anEvent);
}

#endif
