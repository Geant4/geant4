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
// $Id: G4Event.cc,v 1.5 2002-08-13 18:17:53 asaim Exp $
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
 thePrimaryVertex(0),numberOfPrimaryVertex(0),
 HC(0),DC(0),trajectoryContainer(0),eventAborted(false)
{;}

G4Event::G4Event(G4int evID)
:eventID(evID),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),
 HC(0),DC(0),trajectoryContainer(0),eventAborted(false)
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

G4int G4Event::operator==(const G4Event &right) const
{
  return ( eventID == right.eventID );
}

G4int G4Event::operator!=(const G4Event &right) const
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
    for(G4int i=0;i<n_traj;i++)
    { (*trajectoryContainer)[i]->DrawTrajectory(); }
  }

  if(HC)
  {
    G4int n_HC = HC->GetCapacity();
    for(G4int j=0;j<n_HC;j++)
    {
      G4VHitsCollection * VHC = HC->GetHC(j);
      if(VHC) VHC->DrawAllHits();
    }
  }

  if(DC)
  {
    G4int n_DC = DC->GetCapacity();
    for(G4int j=0;j<n_DC;j++)
    {
      G4VDigiCollection * VDC = DC->GetDC(j);
      if(VDC) VDC->DrawAllDigi();
    }
  }
}

