// This code implementation is the intellectual property of
// the RD44 GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: G4TrajectoriesModel.cc,v 1.5 1999-05-26 02:11:42 stanaka Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 
// John Allison  26th August 1998.
// Model which knows how to draw GEANT4 trajectories.

#include "G4TrajectoriesModel.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"

G4TrajectoriesModel::~G4TrajectoriesModel () {}

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
	// GB : default representation should be the most efficent one.
        // Representation with a circle at each step should be an option 
        // choosable from a command.
	// (*TC) [iT] -> DrawTrajectory (50);
	(*TC) [iT] -> DrawTrajectory ();
      }
    }
  }
}
