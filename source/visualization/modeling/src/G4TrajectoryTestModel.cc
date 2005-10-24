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
// $Id: G4TrajectoryTestModel.cc,v 1.1 2005-10-24 11:20:18 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005

#include "G4TrajectoryTestModel.hh"

#include "G4VTrajectory.hh"

G4TrajectoryTestModel::G4TrajectoryTestModel
(const G4String& name,
const G4String& commandPrefix) 
{
  fGlobalTag = name;
  fGlobalDescription =
    "G4TrajectoryTestModel: " + fGlobalTag +
    ", command prefix: " + commandPrefix;

  // Create messenger, passing commandPrefix.
}

G4TrajectoryTestModel::~G4TrajectoryTestModel() {}

void G4TrajectoryTestModel::SetTrajectory
(const G4VTrajectory* trajectory, G4int i_mode)
{
  fpTrajectory = trajectory;
  fI_mode = i_mode;
}

void G4TrajectoryTestModel::DescribeYourselfTo(G4VGraphicsScene& sceneHandler)
{
  // Check pointer to protect against mis-use.
  if (fpTrajectory) {
    G4cout <<
      "G4TrajectoryTestModel::DescribeYourselfTo(G4VGraphicsScene&):"
      "\n  Drawing trajectory with i_mode = " << fI_mode
	   << G4endl;
    fpTrajectory->ShowTrajectory();
    G4cout << G4endl;
   }

  // Zero pointer to protect against future mis-use.
  fpTrajectory = 0;
}
