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
// $Id: G4RTTrackingAction.cc,v 1.3 2001-07-11 10:09:04 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//
//

///////////////////////
//G4RTTrackingAction.cc
///////////////////////


#include "G4RTTrackingAction.hh"
#include "G4RayTrajectory.hh"
#include "G4TrackingManager.hh"
#include "G4ios.hh"


void G4RTTrackingAction :: PreUserTrackingAction(const G4Track* aTrack)
{
  G4RayTrajectory* aTrajectory=new G4RayTrajectory;
  fpTrackingManager->SetTrajectory(aTrajectory);

}







