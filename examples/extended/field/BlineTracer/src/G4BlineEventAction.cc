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
/// \file field/BlineTracer/src/G4BlineEventAction.cc
/// \brief Implementation of the G4BlineEventAction class
//
//
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
  for (size_t i=0; i<fTrajectoryVisAttributes.size(); i++)
    delete fTrajectoryVisAttributes[i];
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::BeginOfEventAction(const G4Event*)
{
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::EndOfEventAction(const G4Event* evt)
{
  G4TrajectoryContainer * trajectoryContainer = evt->GetTrajectoryContainer();
  if(trajectoryContainer) 
  {
    // visualisation
    // -------------
  
    if (fDrawBline || fDrawPoints)
    {
      G4int n_point = (*(evt->GetTrajectoryContainer()))[0]->GetPointEntries();
       
      G4Polyline pPolyline;
      G4Polymarker stepPoints;
      fTrajectoryVisAttributes.push_back(new G4VisAttributes(fDrawColour));
      stepPoints.SetMarkerType(G4Polymarker::circles);
      stepPoints.SetScreenSize(fPointSize);
      stepPoints.SetFillStyle(G4VMarker::filled);
      stepPoints.SetVisAttributes(fTrajectoryVisAttributes.back());

      for(G4int i=0; i<n_point; i++)
      {
        G4ThreeVector pos = ((G4TrajectoryPoint*)
          ((*(evt->GetTrajectoryContainer()))[0]->GetPoint(i)))->GetPosition();
        if (fDrawBline)  pPolyline.push_back( pos);
        if (fDrawPoints)  stepPoints.push_back(pos);
      }

      pPolyline.SetVisAttributes(fTrajectoryVisAttributes.back());

      fTrajectoryPolyline.push_back(pPolyline);
      fTrajectoryPoints.push_back(stepPoints);
    }
  }
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::
DrawFieldLines( G4double, G4double, G4double )
{
  size_t nline = fTrajectoryPolyline.size();
  size_t npoints =fTrajectoryPoints.size();

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
    pVVisManager->Draw(fTrajectoryPolyline[i]);
  for (size_t i=0;i<npoints;i++)
    pVVisManager->Draw(fTrajectoryPoints[i]); 

  // ((G4VisManager*)pVVisManager)->GetCurrentViewer()->DrawView();   
  // ((G4VisManager*)pVVisManager)->GetCurrentViewer()->ShowView();
}

///////////////////////////////////////////////////////////////////////////

void G4BlineEventAction::ResetVectorObjectToBeDrawn()
{
  fTrajectoryVisAttributes.clear();
  fTrajectoryPolyline.clear();
  fTrajectoryPoints.clear();
}
