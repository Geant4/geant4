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
// $Id: G4VTrajectoryModel.hh,v 1.1 2005-10-24 11:20:18 allison Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for trajectory models. Trajectory models are responsible
// for drawing trajectories according to a particular style.
// Class Description - End:

#ifndef G4VTRAJECTORYMODEL_HH
#define G4VTRAJECTORYMODEL_HH

#include "G4VModel.hh"

class G4VTrajectory;

class G4VTrajectoryModel: public G4VModel {

public: // With description

  G4VTrajectoryModel() {}

  virtual ~G4VTrajectoryModel() {}

  virtual void SetTrajectory(const G4VTrajectory*, G4int i_mode = 0) = 0;
  // Set the trajectory with optional i_mode parameter

  virtual void Print() const = 0;
  // Print drawer configuration

protected:

  const G4VTrajectory* fpTrajectory;
  G4int fI_mode;

};

#endif
