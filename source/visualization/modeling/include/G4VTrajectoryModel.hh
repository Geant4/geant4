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
// $Id: G4VTrajectoryModel.hh,v 1.3 2005/12/02 19:16:15 perl Exp $
// GEANT4 tag $Name: geant4-08-00 $
//
// Jane Tinslay, John Allison, Joseph Perl October 2005
//
// Class Description:
// Abstract base class for trajectory models. Trajectory models are responsible
// for drawing individual trajectories according to a particular style.
// Class Description - End:

#ifndef G4VTRAJECTORYMODEL_HH
#define G4VTRAJECTORYMODEL_HH

#include "G4String.hh"
#include <iostream>

class G4VTrajectory;

class G4VTrajectoryModel {

public: // With description

  G4VTrajectoryModel(const G4String& name):fName(name) {}

  virtual ~G4VTrajectoryModel() {}

  virtual void Draw(const G4VTrajectory& trajectory, G4int i_mode = 0) const = 0;
  // Draw the trajectory with optional i_mode parameter

  virtual void Print(std::ostream& ostr) const = 0;
  // Print configuration

  G4String Name() const ;

protected:

  // Data member
  G4String fName;

};

inline G4String G4VTrajectoryModel::Name() const {return fName;}

#endif

