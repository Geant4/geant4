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
//
// $Id: G4RTSteppingAction.cc 74050 2013-09-20 09:38:19Z gcosmo $
//
//
//


#include "G4RTSteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4VisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4TransportationManager.hh"

G4bool G4RTSteppingAction::ignoreTransparency = false;

G4RTSteppingAction::G4RTSteppingAction()
{}

void G4RTSteppingAction::UserSteppingAction(const G4Step* aStep)
{
  // Stop the ray particle if it reaches to the coloured volume with its alpha = 1
  // or coloured volume and the user request to ignore transparency

  G4StepPoint* fpStepPoint = aStep -> GetPostStepPoint();
  G4VPhysicalVolume* postVolume_phys=fpStepPoint->GetPhysicalVolume();

  if(!postVolume_phys) return;   // Reach to out of the world

  const G4LogicalVolume* postVolume_log = postVolume_phys->GetLogicalVolume();
  const G4VisAttributes* postVisAtt = postVolume_log -> GetVisAttributes();
  G4VisManager* visManager = G4VisManager::GetInstance();
  if(visManager) {
    G4VViewer* viewer = visManager->GetCurrentViewer();
    if (viewer) {
      postVisAtt = viewer->GetApplicableVisAttributes(postVisAtt);
    }
  }
  if((!postVisAtt) || (!(postVisAtt->IsVisible()))) return; // Invisible volume, continue tracking

  if((postVisAtt->IsForceDrawingStyle())
    &&(postVisAtt->GetForcedDrawingStyle()==G4VisAttributes::wireframe)) return;
                                          // Wire frame volume, continue tracking

  G4double postAlpha=(postVisAtt->GetColour()).GetAlpha();

  if(postAlpha==1.0 || ignoreTransparency) // Stop stepping
  { 
    G4Track* currentTrack = aStep -> GetTrack();
    currentTrack -> SetTrackStatus(fStopAndKill);
  }
}

void G4RTSteppingAction::SetIgnoreTransparency(G4bool val)
{ ignoreTransparency = val; }
G4bool G4RTSteppingAction::GetIgnoreTransparency()
{ return ignoreTransparency; }
