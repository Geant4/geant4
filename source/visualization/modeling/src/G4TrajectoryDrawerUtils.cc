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
// $Id: G4TrajectoryDrawerUtils.cc,v 1.4 2006-03-24 20:22:43 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005

#include "G4TrajectoryDrawerUtils.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

namespace G4TrajectoryDrawerUtils {

  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints) 
  {
    for (G4int i=0; i<traj.GetPointEntries(); i++) {
      G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
      
      const std::vector<G4ThreeVector>* auxiliaries
	= aTrajectoryPoint->GetAuxiliaryPoints();
      
      if (auxiliaries) {
	for (size_t iAux=0; iAux<auxiliaries->size(); ++iAux) {
	  const G4ThreeVector pos((*auxiliaries)[iAux]);
	  trajectoryLine.push_back(pos);
	  auxiliaryPoints.push_back(pos);
	}
      }
      const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
      trajectoryLine.push_back(pos);
      stepPoints.push_back(pos);
    }    
  }

  void DrawLineAndPoints(const G4VTrajectory& traj, const G4int& i_mode, const G4Colour& colour, const G4bool& visible) {
    // If i_mode>=0, draws a trajectory as a polyline (default is blue for
    // positive, red for negative, green for neutral) and, if i_mode!=0,
    // adds markers - yellow circles for step points and magenta squares
    // for auxiliary points, if any - whose screen size in pixels is     
    // given by std::abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily 
    // visible markers.
    
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (0 == pVVisManager) return;
    
    const G4double markerSize = std::abs(i_mode)/1000;
    G4bool lineRequired (i_mode >= 0);
    G4bool markersRequired (markerSize > 0.);   
    
    // Return if don't need to do anything
    if (!lineRequired && !markersRequired) return;
    
    // Get points to draw
    G4Polyline trajectoryLine;
    G4Polymarker stepPoints;
    G4Polymarker auxiliaryPoints;
    
    GetPoints(traj, trajectoryLine, auxiliaryPoints, stepPoints);
    
    if (lineRequired) {
      G4VisAttributes trajectoryLineAttribs(colour);
      trajectoryLineAttribs.SetVisibility(visible);
      trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
      pVVisManager->Draw(trajectoryLine);
    }
    
    if (markersRequired) {
      auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
      auxiliaryPoints.SetScreenSize(markerSize);
      auxiliaryPoints.SetFillStyle(G4VMarker::filled);
      G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
      auxiliaryPointsAttribs.SetVisibility(visible);
      auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
      pVVisManager->Draw(auxiliaryPoints);
      
      stepPoints.SetMarkerType(G4Polymarker::circles);
      stepPoints.SetScreenSize(markerSize);
      stepPoints.SetFillStyle(G4VMarker::filled);  
      G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
      stepPoints.SetVisAttributes(&stepPointsAttribs);
      pVVisManager->Draw(stepPoints);
    }
    
  }
}
