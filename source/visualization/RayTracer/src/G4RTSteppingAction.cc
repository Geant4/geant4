
#include "G4RTSteppingAction.hh"

#include "G4SteppingManager.hh"
#include "G4VisManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Track.hh"
#include "G4StepStatus.hh"
#include "G4TransportationManager.hh"

G4RTSteppingAction::G4RTSteppingAction()
{
  ignoreTransparency=false;
}

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
  if(visManager) postVisAtt = visManager->GetCurrentViewer()->GetApplicableVisAttributes(postVisAtt);

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

