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
// $Id: G4TrajectoryDrawByCharge.cc,v 1.1 2005-11-02 00:41:13 tinslay Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// Jane Tinslay October 2005

#include "G4TrajectoryDrawByCharge.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"
#include "G4TrajectoryDrawerUtils.hh"
#include "G4VTrajectory.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"

G4TrajectoryDrawByCharge::G4TrajectoryDrawByCharge(const G4String& name)
  :G4VTrajectoryDrawer(name)
  ,fPositive(G4Color::Blue)
  ,fNegative(G4Colour::Red)
  ,fNeutral(G4Colour::Green)
{}

G4TrajectoryDrawByCharge::G4TrajectoryDrawByCharge(const G4Colour& positive,
						   const G4Colour& negative,
						   const G4Colour& neutral)
  :fPositive(positive)
  ,fNegative(negative)
  ,fNeutral(neutral)
{}

G4TrajectoryDrawByCharge::~G4TrajectoryDrawByCharge() {}

void
G4TrajectoryDrawByCharge::Draw(const G4VTrajectory& traj, G4int i_mode) 
{
  // If i_mode>=0, draws a trajectory as a polyline (default is blue for
  // positive, red for negative, green for neutral) and, if i_mode!=0,
  // adds markers - yellow circles for step points and magenta squares
  // for auxiliary points, if any - whose screen size in pixels is     
  // given by std::abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily 
  // visible markers.

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;

  const G4double markerSize = std::abs(i_mode)/1000;
  G4bool lineRequired (i_mode >= 0);
  G4bool markersRequired (markerSize > 0.);   

  // Return if don't need to do anything
  if (!lineRequired && !markersRequired) return;

  // Get points to draw
  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;
  
  G4TrajectoryDrawerUtils::GetPoints(traj, trajectoryLine, 
				     auxiliaryPoints, stepPoints);
  
  if (lineRequired) {
    G4Colour colour;
    const G4double charge = traj.GetCharge();
    if(charge>0.)      colour = fPositive; 
    else if(charge<0.) colour = fNegative; 
    else               colour = fNeutral; 
    
    G4VisAttributes trajectoryLineAttribs(colour);
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    pVVisManager->Draw(trajectoryLine);
  }

  if (markersRequired) {
    auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    auxiliaryPoints.SetScreenSize(markerSize);
    auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    pVVisManager->Draw(auxiliaryPoints);
    
    stepPoints.SetMarkerType(G4Polymarker::circles);
    stepPoints.SetScreenSize(markerSize);
    stepPoints.SetFillStyle(G4VMarker::filled);  
    G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow.
    stepPoints.SetVisAttributes(&stepPointsAttribs);
    pVVisManager->Draw(stepPoints);
  }
  
  return;
}

void
G4TrajectoryDrawByCharge::Print() const
{
  G4cout<<"G4TrajectoryDrawByCharge drawer with name " <<GetName()<<G4endl;
  G4cout<<"Trajectory colour scheme: "<<G4endl;
  G4cout<<"Positive: "<< fPositive <<G4endl;
  G4cout<<"Negative: "<< fNegative <<G4endl;
  G4cout<<"Neutral : "<< fNeutral <<G4endl;
}
