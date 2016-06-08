// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PEvent.cc,v 1.11 1999/12/15 14:51:20 gunter Exp $
// GEANT4 tag $Name: geant4-02-00 $
//

// G4PEvent

#include "G4PEvent.hh"

#include "G4Event.hh"

G4PEvent::G4PEvent()
:eventID(0),
 thePrimaryVertex(NULL),numberOfPrimaryVertex(0),HC(NULL),DC(NULL)
{;}

G4PEvent::G4PEvent(const G4Event *evt)
:eventID(0),
 thePrimaryVertex(NULL),numberOfPrimaryVertex(0),HC(NULL),DC(NULL)
{
  InitPEvent(evt, NULL, NULL);
}

G4PEvent::G4PEvent(const G4Event *evt,
                   HepRef(G4PHCofThisEvent) pHC,
                   HepRef(G4PDCofThisEvent) pDC)
:eventID(0),
 thePrimaryVertex(NULL),numberOfPrimaryVertex(0),HC(NULL),DC(NULL)
{
  InitPEvent(evt, pHC, pDC);
}

void G4PEvent::InitPEvent(const G4Event *evt,
                          HepRef(G4PHCofThisEvent) pHC,
                          HepRef(G4PDCofThisEvent) pDC)
{
  eventID = evt->GetEventID();

  const G4PrimaryVertex* PV = evt->GetPrimaryVertex();
  if(PV) thePrimaryVertex = new(ooThis()) G4PPrimaryVertex(PV);
  numberOfPrimaryVertex = evt->GetNumberOfPrimaryVertex();

  if(pHC!=NULL) HC = pHC;

  if(pDC!=NULL) DC = pDC;

//  const G4TrajectoryContainer* TC = evt->GetTrajectoryContainer();
//  if(TC) trajectoryContainer = new G4PTrajectoryContainer(TC);

}

G4PEvent::~G4PEvent()
{
  if(thePrimaryVertex != NULL)
  {
    G4PPrimaryVertex* aPPV = (HepRef(G4PPrimaryVertex)) thePrimaryVertex;
    HepDelete(aPPV);
  }

  if(HC != NULL)
  {
    G4PHCofThisEvent* aHC = (HepRef(G4PHCofThisEvent)) HC;
    HepDelete(aHC);
  }

  if(DC != NULL)
  {
    G4PDCofThisEvent* aDC = (HepRef(G4PDCofThisEvent)) DC;
    HepDelete(aDC);
  }

//  if(trajectoryContainer)
//    HepDelete(trajectoryContainer);
}

G4Event* G4PEvent::MakeTransientObject()
{
  G4Event* anEvt = new G4Event(eventID);

  if(thePrimaryVertex != NULL)
  {
    G4PrimaryVertex* aPV = thePrimaryVertex->MakeTransientObject();
    for( G4int i = 0; i<numberOfPrimaryVertex; i++ )
    {
      if( aPV != NULL )
      {
        anEvt->AddPrimaryVertex(aPV);
        aPV = aPV->GetNext();
      }
    }
  }

  // Associate persistent Hits Collctions and Digits Collections
  // to transient G4Event ... this requires an introduction of
  // G4VHCofThisEvent and G4VDCofThisEvent in G4Event
  // (will be implemented in next release)
  // if(HC != NULL)
  //   anEvt->SetHCofThisEvent( HC );
  // if(DC != NULL)
  //   anEvt->SetDCofThisEvent( DC );

  return anEvt;
}

void G4PEvent::SetEventID(const G4Event *evt)
{
  eventID = evt->GetEventID();
}

