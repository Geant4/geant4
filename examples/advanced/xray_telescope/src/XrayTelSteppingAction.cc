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
// **********************************************************************
// *                                                                    *
// *                    GEANT 4 xray_telescope advanced example         *
// *                                                                    *
// * MODULE:            XrayTelSteppingAction.cc                        *
// * -------                                                            *
// *                                                                    *
// * Version:           0.4                                             *
// * Date:              06/11/00                                        *
// * Author:            R Nartallo                                      *
// * Organisation:      ESA/ESTEC, Noordwijk, THe Netherlands           *
// *                                                                    *
// **********************************************************************
// 
// CHANGE HISTORY
// --------------
//
// 05.12.2001 R. Nartallo
// - Return condition for entering detector (cf LowEnTest)
// - Remove line to kill track if too many steps
//
// 07.11.2001 M.G. Pia
// - Modified the analysis management
// - Small design iteration
//
// 30.11.2000 R. Nartallo
// - Add pre-processor directives to compile without analysis option
//
// 16.11.2000 A. Pfeiffer
// - Implementation of analysis manager call
//
// 06.11.2000 R.Nartallo
// - First implementation of xray_telescope Physics list
// - Based on Chandra and XMM models
// 
// **********************************************************************

#include "XrayTelSteppingAction.hh"
#include "globals.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "XrayTelEventAction.hh"
#include "XrayTelRunAction.hh"
#include "XrayTelAnalysis.hh"

#include "Randomize.hh"

XrayTelSteppingAction::XrayTelSteppingAction()
{ }


XrayTelSteppingAction::~XrayTelSteppingAction()
{ }


void XrayTelSteppingAction::UserSteppingAction(const G4Step* step)
{
  G4bool entering = false;
  G4Track* track = step->GetTrack();

  G4String volName; 
  if (track->GetVolume()) volName =  track->GetVolume()->GetName(); 
  G4String nextVolName;
  if (track->GetNextVolume()) nextVolName =  track->GetNextVolume()->GetName();

  // Entering Detector
  if (volName != "Detector_P" && nextVolName == "Detector_P") 
    {
      entering = true;

      // Notify the corresponding UserAction that the event must be visualised
      G4RunManager* runManager = G4RunManager::GetRunManager();
      // Casting is safe here: one knows the RTTI of the UserActions 
      // in the current application (otherwise one could dynamic_cast)
      XrayTelEventAction* eventAction = (XrayTelEventAction*) runManager->GetUserEventAction();
      eventAction->Update();

      // Notify the corresponding UserAction to update the run counters
      XrayTelRunAction* runAction = (XrayTelRunAction*) runManager->GetUserRunAction();
      runAction->Update(track->GetKineticEnergy());
    }

  // Do the analysis related to this step
  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
  analysis->analyseStepping(*track,entering);
}




