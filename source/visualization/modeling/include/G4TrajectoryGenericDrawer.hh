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
// $Id: G4TrajectoryGenericDrawer.hh,v 1.1 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description:
// Trajectory model which colours a trajectory according to  
// charge. Guts taken from G4VTrajectory::DrawTrajectory method.
// Class Description - End:

#ifndef G4TRAJECTORYGENERICDRAWER_HH
#define G4TRAJECTORYGENERICDRAWER_HH

#include "G4Colour.hh"
#include "G4VTrajectoryModel.hh"

class G4TrajectoryGenericDrawer : public G4VTrajectoryModel {

public: // With description

  G4TrajectoryGenericDrawer(const G4String& name = "Unspecified", G4VisTrajContext* context=0);

  virtual ~G4TrajectoryGenericDrawer();

  virtual void Draw(const G4VTrajectory& trajectory, const G4int& i_mode = 0, 
		    const G4bool& visible = true) const;
  // Draw the trajectory with optional i_mode parameter

  virtual void Print(std::ostream& ostr) const;
  // Print configuration

private:
  
};

#endif
