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
// $Id: G4TrajectoryDrawerUtils.hh,v 1.2 2005/11/21 05:44:44 tinslay Exp $
// GEANT4 tag $Name: geant4-08-00 $
// 
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
// Class Description
// Simple trajectory model drawing utility functions
// Class Description - End:

#ifndef G4TRAJECTORYDRAWERUTILS
#define G4TRAJECTORYDRAWERUTILS

class G4Polyline;
class G4Polymarker;
class G4VTrajectory;

namespace G4TrajectoryDrawerUtils {

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints);

}

#endif
