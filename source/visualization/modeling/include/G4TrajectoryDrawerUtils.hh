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
// $Id: G4TrajectoryDrawerUtils.hh,v 1.5 2006-05-02 20:47:40 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description
// Simple trajectory model drawing utility functions
// Class Description - End:

#ifndef G4TRAJECTORYDRAWERUTILS
#define G4TRAJECTORYDRAWERUTILS

#include "globals.hh"

class G4Colour;
class G4Polyline;
class G4Polymarker;
class G4VisTrajContext;
class G4VTrajectory;

namespace G4TrajectoryDrawerUtils {

  namespace {
    G4bool warnedAboutIMode = false;
  }

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints);

  // Draw trajectory line and points
  void DrawLineAndPoints(const G4VTrajectory& traj, const G4int& i_mode, 
			 const G4Colour& cfg, const G4bool& visible);

  // Draw trajectory line and points using G4VisTrajContext object information
  void DrawLineAndPoints(const G4VTrajectory& traj, const G4VisTrajContext&, const G4int& i_mode=0);

}

#endif
