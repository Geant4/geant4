// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoriesModel.cc,v 1.1 1999-01-07 16:15:38 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 trajectories.

#include "G4TrajectoriesModel.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

G4TrajectoriesModel::G4TrajectoriesModel () {
  fGlobalTag = "G4TrajectoriesModel for all trajectories.";
  fGlobalDescription = fGlobalTag;
}

void G4TrajectoriesModel::DescribeYourselfTo (G4VGraphicsScene& scene) {
  G4RunManager* runManager = G4RunManager::GetRunManager ();
  const G4Event* event = runManager -> GetCurrentEvent ();
  if (event) {
    G4TrajectoryContainer* TC = event -> GetTrajectoryContainer ();
    if (TC) {
      G4int nT = TC -> entries ();
      for (int iT = 0; iT < nT; iT++) {
	(*TC) [iT] -> DrawTrajectory (50);
      }
    }
  }
}
