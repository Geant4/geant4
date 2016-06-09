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
// $Id: G4TrajectoryDrawerUtils.cc,v 1.6 2006/06/29 21:33:12 gunter Exp $
// GEANT4 tag $Name: geant4-08-01 $
//
// Jane Tinslay, John Allison, Joseph Perl November 2005
//
#include "G4TrajectoryDrawerUtils.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4VisAttributes.hh"
#include "G4VisTrajContext.hh"
#include "G4VVisManager.hh"

namespace G4TrajectoryDrawerUtils {

  
  void GetPoints(const G4VTrajectory& traj, G4Polyline& trajectoryLine,
		 G4Polymarker& auxiliaryPoints, G4Polymarker& stepPoints) 
  {
    for (G4int i=0; i<traj.GetPointEntries(); i++) {
      G4VTrajectoryPoint* aTrajectoryPoint = traj.GetPoint(i);
      
      const std::vector<G4ThreeVector>* auxiliaries
	= aTrajectoryPoint->GetAuxiliaryPoints();
      
      if (0 != auxiliaries) {
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
  
  void DrawLineAndPoints(const G4VTrajectory& traj, const G4VisTrajContext& context, const G4int& i_mode) 
  {
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if (0 == pVVisManager) return;
    
    // Extra copy while i_mode is still around
    G4VisTrajContext myContext(context);
    
    if (i_mode != 0) {
      const G4double markerSize = std::abs(i_mode)/1000;
      G4bool lineRequired (i_mode >= 0);
      G4bool markersRequired (markerSize > 0.);        
      
      myContext.SetDrawLine(lineRequired);
      myContext.SetDrawAuxPts(markersRequired);
      myContext.SetDrawStepPts(markersRequired);

      myContext.SetAuxPtsSize(markerSize);
      myContext.SetStepPtsSize(markerSize);

      if (!warnedAboutIMode) {
	G4cout<<"Trajectory drawing configuration will be based on imode value of "<<i_mode<<G4endl;
	warnedAboutIMode = true;
      }
    }

    // Return if don't need to do anything
    if (!myContext.GetDrawLine() && !myContext.GetDrawAuxPts() && !myContext.GetDrawStepPts()) return;
    
    // Get points to draw
    G4Polyline trajectoryLine;
    G4Polymarker stepPoints;
    G4Polymarker auxiliaryPoints;
    
    GetPoints(traj, trajectoryLine, auxiliaryPoints, stepPoints);
    
    if (myContext.GetDrawLine()) {
      G4VisAttributes trajectoryLineAttribs(myContext.GetLineColour());
      trajectoryLineAttribs.SetVisibility(myContext.GetLineVisible());
      trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
      
      pVVisManager->Draw(trajectoryLine);
    }
    
    if (myContext.GetDrawAuxPts() && (auxiliaryPoints.size() > 0)) {
      auxiliaryPoints.SetMarkerType(myContext.GetAuxPtsType());
      auxiliaryPoints.SetScreenSize(myContext.GetAuxPtsSize());
      auxiliaryPoints.SetFillStyle(myContext.GetAuxPtsFillStyle());
      
      G4VisAttributes auxiliaryPointsAttribs(myContext.GetAuxPtsColour());  
      auxiliaryPointsAttribs.SetVisibility(myContext.GetAuxPtsVisible());
      auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
      
      pVVisManager->Draw(auxiliaryPoints);
    }
    
    if (myContext.GetDrawStepPts() && (stepPoints.size() > 0)) {
      stepPoints.SetMarkerType(myContext.GetStepPtsType());
      stepPoints.SetScreenSize(myContext.GetStepPtsSize());
      stepPoints.SetFillStyle(myContext.GetStepPtsFillStyle());
      
      G4VisAttributes stepPointsAttribs(myContext.GetStepPtsColour()); 
      stepPointsAttribs.SetVisibility(myContext.GetStepPtsVisible());
      stepPoints.SetVisAttributes(&stepPointsAttribs);
      
      pVVisManager->Draw(stepPoints);
    }
  }
}
