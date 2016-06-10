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
#include "G4Threading.hh"

#include "XrayTelAnalysis.hh"

XrayTelSteppingAction::XrayTelSteppingAction() 
{;}


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

  XrayTelAnalysis* analysis = XrayTelAnalysis::getInstance();
 
  // Entering Detector
  if (volName != "Detector_P" && nextVolName == "Detector_P") 
    {
      entering = true;
      analysis->Update(track->GetKineticEnergy(),G4Threading::G4GetThreadId());
    }

  // Do the analysis related to this step
 
  analysis->analyseStepping(*track,entering);

}




