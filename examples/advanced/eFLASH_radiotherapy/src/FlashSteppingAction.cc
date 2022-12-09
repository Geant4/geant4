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
/// \file FlashSteppingAction.cc
/// \brief Implementation of the FlashSteppingAction class

#include "FlashSteppingAction.hh"
#include "FlashDetectorConstruction.hh"
#include "FlashEventAction.hh"
#include "FlashRunAction.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4LogicalVolume.hh"
#include "G4RunManager.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include <sstream>
#include <string>

FlashSteppingAction::FlashSteppingAction(FlashEventAction *)
    : G4UserSteppingAction()
     
{


}

FlashSteppingAction::~FlashSteppingAction() {
}

void FlashSteppingAction::UserSteppingAction(const G4Step *aStep)
{

  G4int eventid =
      G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  G4StepPoint *postStep = aStep->GetPostStepPoint();
  G4StepPoint *preStep = aStep->GetPreStepPoint();
  // G4int trackID = aStep->GetTrack()->GetTrackID();

  //===========================INCIDENT
  //PARTICLES===============================================

  //==============================PRIMARIES ENERGY
  //SPECTRUM=======================================
  if (postStep->GetStepStatus() == fGeomBoundary)
    {

    G4String volumeName =
        postStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();
    G4String prevolumeName =
        preStep->GetPhysicalVolume()->GetLogicalVolume()->GetName();

 G4int parentid = aStep->GetTrack()->GetParentID();
 //G4Track* theTrack = aStep->GetTrack();
 G4double kineticEnergy = aStep->GetTrack()->GetKineticEnergy();
 G4double pos_x = aStep->GetTrack()->GetPosition().x();
 G4double pos_y = aStep->GetTrack()->GetPosition().y();
 G4double pos_z = aStep->GetTrack()->GetPosition().z();
 G4double cos_x = aStep->GetTrack()->GetMomentum().x();
 G4double cos_y = aStep->GetTrack()->GetMomentum().y();
 G4double cos_z = aStep->GetTrack()->GetMomentum().z();
 
 G4double momentum = std::sqrt(cos_x*cos_x+cos_y*cos_y+cos_z*cos_z);
 


 if (aStep->GetTrack()->GetDefinition()==G4Electron::ElectronDefinition() && prevolumeName == "logicTreatmentRoom" && volumeName == "phantomLog")

     	{

	  std::ofstream WriteDataIn("PhantomEntrance.txt", std::ios::app);
	  WriteDataIn	 

	    << parentid << "\t"<< eventid << "\t" << kineticEnergy << "\t" <<pos_x<<"\t"<<pos_y<<"\t"<<pos_z<<"\t"<<momentum<<"\t"<<cos_x<<"\t"<<cos_y<<"\t"<<cos_z<<"\t"<<G4endl;
       
	
	 }
	   
    }
}
