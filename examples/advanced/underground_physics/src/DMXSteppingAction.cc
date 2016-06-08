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
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: a.s.howard@ic.ac.uk
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// History:
// 21 Feb 2002 AH: Added Analysis
//
// SteppingAction program
// --------------------------------------------------------------

#include "DMXSteppingAction.hh"
#include "DMXSteppingActionMessenger.hh"

#include "DMXEventAction.hh"

#ifdef G4ANALYSIS_USE
#include "DMXAnalysisManager.hh"
#endif

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4Polyline.hh" 
#include "G4VisAttributes.hh"
#include "globals.hh"
#include "G4ios.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXSteppingAction::DMXSteppingAction(DMXEventAction* eventAction)
  : evtAction(eventAction)  {

  steppingMessenger = new DMXSteppingActionMessenger(this);

  // defaults for messenger
  colourNeutronFlag      = "magenta";
  colourGammaFlag        = "cyan";
  colourOpticalFlag      = "white";
  colourChargedPlusFlag  = "red";
  colourChargedMinusFlag = "blue";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXSteppingAction::~DMXSteppingAction()
{

  delete steppingMessenger;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXSteppingAction::UserSteppingAction(const G4Step* fStep)
{

  // removed 28/11/01 - unnecessary unless program "freezes"
  // kill track if too many steps
  // NB: This is set to DBL_MAX - therefore may cause program to "hang"
  //  G4int MaxNoSteps = DBL_MAX;
  //  G4int StepNo = fStep->GetTrack()->GetCurrentStepNumber();
  //  if(StepNo >= MaxNoSteps) fStep->GetTrack()->SetTrackStatus(fStopAndKill);

#ifdef G4ANALYSIS_USE 
  G4int StepNo = fStep->GetTrack()->GetCurrentStepNumber();
  if(StepNo == 1) 
    { 
      G4double partEnergy = fStep->GetPreStepPoint()->GetKineticEnergy();
      G4ParticleDefinition* particleType = fStep->GetTrack()->GetDefinition();
      G4String particleName = particleType->GetParticleName();
      DMXAnalysisManager* analysis =  DMXAnalysisManager::getInstance();
      analysis->analyseParticleSource(partEnergy, particleName);
    }
#endif

  // check what is to be drawn from EventAction/EventActionMessenger
  G4String drawColsFlag = evtAction->GetDrawColsFlag();
  G4String drawTrksFlag = evtAction->GetDrawTrksFlag();

  // draw by step (here) instead of by event (event action)
  if (drawColsFlag=="custom" && drawTrksFlag!="none") {

    // check that VisManager exists
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager) {

      // particle colour in a string
      G4String name = fStep->GetTrack()->GetDefinition()->GetParticleName();
      G4String strColour;
      if(name=="neutron") {
	if(drawTrksFlag=="charged") return;
	strColour = colourNeutronFlag;
      } else if (name=="gamma") {
	if(drawTrksFlag=="charged") return;
	strColour = colourGammaFlag;
      } else if (name=="opticalphoton") {
	if(drawTrksFlag!="all") return;
	strColour = colourOpticalFlag;
      }
      else if (name=="alpha" || name=="e+")
	strColour = colourChargedPlusFlag;
      else
	strColour = colourChargedMinusFlag;

      // convert string to G4Colour
      G4Colour colour;
      if     (strColour=="white")    colour=G4Colour(1.0, 1.0, 1.0);
      else if(strColour=="grey" )    colour=G4Colour(0.5, 0.5, 0.5);
      else if(strColour=="lgrey")    colour=G4Colour(.75, .75, .75);
      else if(strColour=="black")    colour=G4Colour(0.0, 0.0, 0.0);
      else if(strColour=="red")      colour=G4Colour(1.0, 0.0, 0.0);
      else if(strColour=="green")    colour=G4Colour(0.0, 1.0, 0.0);
      else if(strColour=="blue")     colour=G4Colour(0.0, 0.0, 1.0);
      else if(strColour=="cyan")     colour=G4Colour(0.0, 1.0, 1.0);
      else if(strColour=="magenta")  colour=G4Colour(1.0, 0.0, 1.0);
      else if(strColour=="yellow")   colour=G4Colour(1.0, 1.0, 0.0);
      else if(strColour=="lgreen")   colour=G4Colour(0.0, .75, 0.0);
      else if(strColour=="lblue")    colour=G4Colour(0.0, 0.0, .75);
      else                           colour=G4Colour(1.0, 1.0, 1.0);

      // create line with colour
      G4VisAttributes attribs(colour);
      G4Polyline polyline;
      polyline.SetVisAttributes(attribs);

      // draw line
      G4Point3D start(fStep->GetPreStepPoint()->GetPosition());
      G4Point3D end(fStep->GetPostStepPoint()->GetPosition());
      polyline.push_back(start);
      polyline.push_back(end);
      pVVisManager->Draw(polyline);
    }
    
  }
  
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

