// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4Event.cc,v 1.2 1999-12-15 14:49:40 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

// G4Event

#include "G4Event.hh"
#include "G4VVisManager.hh"
//#include "G4HCofThisEvent.hh"
//#include "G4DCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4VDigiCollection.hh"
#include "G4ios.hh"

G4Allocator<G4Event> anEventAllocator;

G4Event::G4Event()
:eventID(0),
 thePrimaryVertex(NULL),numberOfPrimaryVertex(0),
 DC(NULL),HC(NULL),trajectoryContainer(NULL)
{;}

G4Event::G4Event(G4int evID)
:eventID(evID),
 thePrimaryVertex(NULL),numberOfPrimaryVertex(0),
 DC(NULL),HC(NULL),trajectoryContainer(NULL)
{;}

G4Event::~G4Event()
{ 
  if(thePrimaryVertex) delete thePrimaryVertex;
  if(HC) delete HC;
  if(DC) delete DC;
  if(trajectoryContainer)
  {
    trajectoryContainer->clearAndDestroy();
    delete trajectoryContainer;
  }
}

int G4Event::operator==(const G4Event &right) const
{
  return ( eventID == right.eventID );
}

int G4Event::operator!=(const G4Event &right) const
{
  return ( eventID != right.eventID );
}

void G4Event::Print() const
{
  G4cout << "G4Event " << eventID << G4endl;
}

void G4Event::Draw() const
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(!pVVisManager) return;

  if(trajectoryContainer)
  {
    G4int n_traj = trajectoryContainer->entries();
    for(int i=0;i<n_traj;i++)
    { (*trajectoryContainer)[i]->DrawTrajectory(); }
  }

  if(HC)
  {
    G4int n_HC = HC->GetCapacity();
    for(int j=0;j<n_HC;j++)
    {
      G4VHitsCollection * VHC = HC->GetHC(j);
      if(VHC) VHC->DrawAllHits();
    }
  }

  if(DC)
  {
    G4int n_DC = DC->GetCapacity();
    for(int j=0;j<n_DC;j++)
    {
      G4VDigiCollection * VDC = DC->GetDC(j);
      if(VDC) VDC->DrawAllDigi();
    }
  }
}

