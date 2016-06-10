//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: G4TrajectoryDrawerUtils.hh 66373 2012-12-18 09:41:34Z gcosmo $
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

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints);

  // Draw trajectory line and points using G4VisTrajContext object information
  void DrawLineAndPoints(const G4VTrajectory& traj, const G4VisTrajContext&);

}

#endif
