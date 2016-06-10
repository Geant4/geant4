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
// $Id: G4Event.cc 84194 2014-10-10 14:29:58Z gcosmo $
//

// G4Event

#include "G4Event.hh"
#include "G4VVisManager.hh"
//#include "G4HCofThisEvent.hh"
//#include "G4DCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4VDigiCollection.hh"
#include "G4ios.hh"

G4ThreadLocal G4Allocator<G4Event> *anEventAllocator = 0;

G4Event::G4Event()
:eventID(0),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),
 HC(0),DC(0),trajectoryContainer(0),eventAborted(false),userInfo(0),
 randomNumberStatus(0),validRandomNumberStatus(false),
 randomNumberStatusForProcessing(0),validRandomNumberStatusForProcessing(false),
 keepTheEvent(false),grips(0)
{
}

G4Event::G4Event(G4int evID)
:eventID(evID),
 thePrimaryVertex(0),numberOfPrimaryVertex(0),
 HC(0),DC(0),trajectoryContainer(0),eventAborted(false),userInfo(0),
 randomNumberStatus(0),validRandomNumberStatus(false),
 randomNumberStatusForProcessing(0),validRandomNumberStatusForProcessing(false),
 keepTheEvent(false),grips(0)
{
}

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
  if(userInfo) delete userInfo;
  if(validRandomNumberStatus) delete randomNumberStatus;
  if(validRandomNumberStatusForProcessing) delete randomNumberStatusForProcessing;
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
