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
// $Id: G4TrajectoryTestModel.hh,v 1.2 2005-10-24 14:03:36 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Test model for drawing trajectories.
// Class Description - End:

#ifndef G4TRAJECTORYTESTMODEL_HH
#define G4TRAJECTORYTESTMODEL_HH

#include "globals.hh"
#include "G4String.hh"
#include "G4VTrajectoryModel.hh"

class G4VTrajectory;

class G4TrajectoryTestModel: public G4VTrajectoryModel {

public: // With description

  G4TrajectoryTestModel
  (const G4String& name,
   const G4String& commandPrefix = "/");

  ~G4TrajectoryTestModel();

  void SetTrajectory(const G4VTrajectory*, G4int i_mode = 0);
  // Set the trajectory with optional i_mode parameter

  void DescribeYourselfTo(G4VGraphicsScene&);

  void Print() const {}
  // Print drawer configuration

};

#endif
