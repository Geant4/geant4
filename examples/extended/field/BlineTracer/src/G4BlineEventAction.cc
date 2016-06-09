//
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
//
// $Id: G4BlineEventAction.cc,v 1.1 2003/11/25 09:29:47 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//
// --------------------------------------------------------------------
//
// G4BlineEventAction implementation
//
// --------------------------------------------------------------------
// Author: Laurent Desorgher (desorgher@phim.unibe.ch)
//         Created - 2003-10-06
// --------------------------------------------------------------------

#include "G4BlineEventAction.hh"
#include "G4BlineTracer.hh"

#include "G4Event.hh"
#include "G4Trajectory.hh"
#include "G4EventManager.hh"
#include "G4VisManager.hh"
#include "G4UImanager.hh"
#include "G4Polyline.hh"
#include "G4Polymarker.hh"

///////////////////////////////////////////////////////////////////////////

G4BlineEventAction::G4BlineEventAction(G4BlineTracer* aBlineTool)
{
  fBlineTool=aBlineTool;
}

///////////////////////////////////////////////////////////////////////////

G4BlineEventAction::~G4BlineEventAction()
{
  for (size_t i=0; i<TrajectoryVisAttributes.size(); i++)
    delete TrajectoryVisAttributes[i];
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::BeginOfEventAction(const G4Event*)
{
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::EndOfEventAction(const G4Event* evt)
{
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if(trajectoryContainer) 
  {
    n_trajectories = trajectoryContainer->entries(); 
   
    // visualisation
    // -------------
  
    if (DrawBline || DrawPoints)
    {
      G4int n_point = (*(evt->GetTrajectoryContainer()))[0]->GetPointEntries();
       
      G4Polyline pPolyline;
      G4Polymarker stepPoints;
      TrajectoryVisAttributes.push_back(new G4VisAttributes(DrawColour));
      stepPoints.SetMarkerType(G4Polymarker::circles);
      stepPoints.SetScreenSize(PointSize);
      stepPoints.SetFillStyle(G4VMarker::filled);
      stepPoints.SetVisAttributes(TrajectoryVisAttributes.back());

      for(G4int i=0; i<n_point; i++)
      {
        G4ThreeVector pos = ((G4TrajectoryPoint*)
          ((*(evt->GetTrajectoryContainer()))[0]->GetPoint(i)))->GetPosition();
        if (DrawBline)  pPolyline.push_back( pos);
        if (DrawPoints)  stepPoints.push_back(pos);
      }

      pPolyline.SetVisAttributes(TrajectoryVisAttributes.back());

      TrajectoryPolyline.push_back(pPolyline);
      TrajectoryPoints.push_back(stepPoints);
    }
  }
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::
DrawFieldLines( G4double, G4double, G4double )
{
  size_t nline = TrajectoryPolyline.size();
  size_t npoints =TrajectoryPoints.size();

  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) 
  {
    G4Exception("G4BlineEventAction::DrawFieldLines()",
        "NullPointer", JustWarning,
        "Missing visualisation driver for visualising magnetic field lines!");
    return;
  }
     
  if (nline ==0)
  {
    G4cout << "WARNING - G4BlineEventAction::DrawFieldLines()" << G4endl
           << "          There is nothing to visualise !" << G4endl;
    return;
  }
  ((G4VisManager*)pVVisManager)->GetCurrentSceneHandler()-> ClearStore ();
  G4UImanager::GetUIpointer () -> ApplyCommand ("/vis/drawVolume");

  for (size_t i=0;i<nline;i++)
    pVVisManager->Draw(TrajectoryPolyline[i]);
  for (size_t i=0;i<npoints;i++)
    pVVisManager->Draw(TrajectoryPoints[i]); 

  // ((G4VisManager*)pVVisManager)->GetCurrentViewer()->DrawView();   
  // ((G4VisManager*)pVVisManager)->GetCurrentViewer()->ShowView();
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::ResetVectorObjectToBeDrawn()
{
  TrajectoryVisAttributes.clear();
  TrajectoryPolyline.clear();
  TrajectoryPoints.clear();
}
