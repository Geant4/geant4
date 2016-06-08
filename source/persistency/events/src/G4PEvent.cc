//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: G4PEvent.cc,v 1.13 2001/07/11 10:02:16 gunter Exp $
// GEANT4 tag $Name: geant4-04-00 $
//

// G4PEvent

#include "G4PEvent.hh"

#include "G4Event.hh"

G4PEvent::G4PEvent()
:eventID(0),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),HC(0),DC(0)
{;}

G4PEvent::G4PEvent(const G4Event *evt)
:eventID(0),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),HC(0),DC(0)
{
  InitPEvent(evt, 0, 0);
}

G4PEvent::G4PEvent(const G4Event *evt,
                   HepRef(G4PHCofThisEvent) pHC,
                   HepRef(G4PDCofThisEvent) pDC)
:eventID(0),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),HC(0),DC(0)
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

  if(pHC!=0) HC = pHC;

  if(pDC!=0) DC = pDC;

//  const G4TrajectoryContainer* TC = evt->GetTrajectoryContainer();
//  if(TC) trajectoryContainer = new G4PTrajectoryContainer(TC);

}

G4PEvent::~G4PEvent()
{
  if(thePrimaryVertex != 0)
  {
    G4PPrimaryVertex* aPPV = (HepRef(G4PPrimaryVertex)) thePrimaryVertex;
    HepDelete(aPPV);
  }

  if(HC != 0)
  {
    G4PHCofThisEvent* aHC = (HepRef(G4PHCofThisEvent)) HC;
    HepDelete(aHC);
  }

  if(DC != 0)
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

  if(thePrimaryVertex != 0)
  {
    G4PrimaryVertex* aPV = thePrimaryVertex->MakeTransientObject();
    for( G4int i = 0; i<numberOfPrimaryVertex; i++ )
    {
      if( aPV != 0 )
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
  // if(HC != 0)
  //   anEvt->SetHCofThisEvent( HC );
  // if(DC != 0)
  //   anEvt->SetDCofThisEvent( DC );

  return anEvt;
}

void G4PEvent::SetEventID(const G4Event *evt)
{
  eventID = evt->GetEventID();
}

