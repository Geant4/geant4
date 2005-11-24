//
// File name:     RadmonApplicationEventTracks.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventTracks.cc,v 1.1 2005-11-24 02:34:21 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonApplicationEventTracks.hh"
#include "G4Event.hh"
#include "G4VVisManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4VTrajectory.hh"

void                                            RadmonApplicationEventTracks :: OnEndOfEvent(const G4Event * event)
{
 if (!G4VVisManager::GetConcreteInstance())
  return;

 G4TrajectoryContainer * trajectoryContainer(event->GetTrajectoryContainer());
 
 if (!trajectoryContainer)
  return;
  
 G4int nTrajectories(trajectoryContainer->entries());
 while (nTrajectories)
 {
  nTrajectories--;
  (*trajectoryContainer)[nTrajectories]->DrawTrajectory(0);
 }
}
