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
// File name:     RadmonApplicationEventTracks.cc
// Creation date: Nov 2005
// Main author:   Riccardo Capra <capra@ge.infn.it>
//
// Id:            $Id: RadmonApplicationEventTracks.cc,v 1.2 2006-06-28 13:45:50 gunter Exp $
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
