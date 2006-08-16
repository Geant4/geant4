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
// $Id: ExN03TrackingAction.cc,v 1.5 2006-08-16 15:44:06 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $

#include "ExN03TrackingAction.hh"
//#include "G4SmoothTrajectory.hh"
#include "G4RichTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4IdentityTrajectoryFilter.hh"
#include "G4TransportationManager.hh"
#include "G4PropagatorInField.hh"

void ExN03TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  /* Instate the following code for advanced trajectory displaying
  // Require rich trajectory...
  fpTrackingManager->SetTrajectory(new G4RichTrajectory(aTrack));

  // Activate storing of auxiliary points for smoother trajectory...
  static G4IdentityTrajectoryFilter curvedFilter;
  G4TransportationManager::GetTransportationManager()->
    GetPropagatorInField()->SetTrajectoryFilter(&curvedFilter);
  */
}
