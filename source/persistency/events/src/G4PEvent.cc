// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4PEvent.cc,v 1.1 1999/01/07 16:10:47 gunter Exp $
// GEANT4 tag $Name: geant4-00-01 $
//

// G4PEvent

#include "G4PEvent.hh"
#include "G4PPrimaryVertex.hh"

G4PEvent::G4PEvent()
{
  eventID =0;
  // initialization of thePrimaryVertex is not necessary
  // since Objectivity/DB will do it for you...
    HepRef(G4PPrimaryVertex) thePrimaryVertex = NULL;
  numberOfPrimaryVertex = 0;
//  HC = NULL;
//  trajectoryContainer = NULL;
}

G4PEvent::G4PEvent(const G4Event *evt)
{
  eventID = evt->GetEventID();

  const G4PrimaryVertex* PV = evt->GetPrimaryVertex();
  if(PV)
  { HepRef(G4PPrimaryVertex) thePrimaryVertex = new G4PPrimaryVertex(PV); }
  else
  {  HepRef(G4PPrimaryVertex) thePrimaryVertex = NULL; }

//  const G4TrajectoryContainer* TC = evt->GetTrajectoryContainer();
//  if(TC)
//  { trajectoryContainer = new G4PTrajectoryContainer(TC); }
//  else
//  { trajectoryContainer = NULL; }
}

G4PEvent::~G4PEvent()
{
//  if(thePrimaryVertex) delete thePrimaryVertex;
//  if(trajectoryContainer) delete trajectoryContainer;
}

