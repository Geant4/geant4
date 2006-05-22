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
// $Id: G4TrajectoriesModel.cc,v 1.17 2006-05-22 07:26:24 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 trajectories.

#include "G4TrajectoriesModel.hh"

#include "G4VGraphicsScene.hh"
#include "G4EventManager.hh"
#include "G4Event.hh"

G4TrajectoriesModel::~G4TrajectoriesModel () {}

G4TrajectoriesModel::G4TrajectoriesModel (G4int drawingMode):
fDrawingMode(drawingMode) {
  fGlobalTag = "G4TrajectoriesModel for all trajectories.";
  fGlobalDescription = fGlobalTag;
}

void G4TrajectoriesModel::DescribeYourselfTo (G4VGraphicsScene& sceneHandler) {
  const G4Event* event = 
    G4EventManager::GetEventManager()->GetConstCurrentEvent();
  if (event) {
    G4TrajectoryContainer* TC = event -> GetTrajectoryContainer ();
    if (TC) {
      for (G4int iT = 0; iT < TC->entries(); iT++) {
	sceneHandler.AddCompound (*((*TC) [iT]));
      }
    }
  }
}
