//
// File name:     RadmonPhysicsSteppingAction.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonPhysicsSteppingAction.cc,v 1.1 2005-12-06 19:36:23 capra Exp $
// Tag:           $Name: not supported by cvs2svn $
//

// Include files
#include "RadmonPhysicsSteppingAction.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

void                                            RadmonPhysicsSteppingAction :: OnUserStepping(const G4Step * step)
{
 G4Track * track(step->GetTrack());
 G4ParticleDefinition * particle(track->GetDefinition());
 G4ProcessManager * manager(particle->GetProcessManager());
 
 if (manager->GetProcessListLength()==0)
  track->SetTrackStatus(fStopAndKill);
}

