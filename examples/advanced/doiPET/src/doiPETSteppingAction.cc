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

//GEANT4 - Depth-of-Interaction enabled Positron emission tomography (PET) advanced example 

//Authors and contributors

// Author list to be updated, with names of co-authors and contributors from National Institute of Radiological Sciences (NIRS)

// Abdella M. Ahmed (1, 2), Andrew Chacon (1, 2), Harley Rutherford (1, 2),
// Hideaki Tashima (3), Go Akamatsu (3), Akram Mohammadi (3), Eiji Yoshida (3), Taiga Yamaya (3)
// Susanna Guatelli (2), and Mitra Safavi-Naeini (1, 2)

// (1) Australian Nuclear Science and Technology Organisation, Australia
// (2) University of Wollongong, Australia
// (3) National Institute of Radiological Sciences, Japan


#include "doiPETSteppingAction.hh"
#include "doiPETAnalysis.hh"
#include "doiPETRun.hh"
#include "doiPETEventAction.hh"

#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4SteppingManager.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <limits>

////////// Constructor //////////////////////////////////////////////
/*doiPETSteppingAction::doiPETSteppingAction()
{;}*/
doiPETSteppingAction::doiPETSteppingAction()
	: G4UserSteppingAction()
{ }

///////// Destructor ////////////////////////////////////////////////
doiPETSteppingAction::~doiPETSteppingAction()
{;}

///////// UserSteppingAction ///////////////////////////////////////
void doiPETSteppingAction::UserSteppingAction(const G4Step* aStep)
{
	doiPETRun* run = static_cast<doiPETRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());

	G4Track* track = aStep->GetTrack();
	G4String volumeName = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
	G4double edep = aStep->GetTotalEnergyDeposit();
	G4String processName =aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	G4String particleName = track->GetDynamicParticle()->GetDefinition()->GetParticleName();
	//G4int eventID = G4RunManager::GetRunManager()->GetCurrentEvent() -> GetEventID();	
	G4ThreeVector pos = track->GetPosition();

	//G4StepPoint* point1 = aStep->GetPreStepPoint();
	//G4ThreeVector posIntCrystal = point1->GetPosition();

	//G4StepPoint* p2 = aStep->GetPostStepPoint();
	G4StepPoint* p1 = aStep->GetPreStepPoint();
	G4ThreeVector coord1 = p1->GetPosition();

	//The following is to get the local position of the with respect to the crystal volume
	const G4AffineTransform transformation =  p1->GetTouchable()-> GetHistory()->GetTopTransform();
	G4ThreeVector localPosition = transformation.TransformPoint(coord1);

	G4int blockID;
	G4int crystalID;

	G4int scatterIndex = 0;//For checking

	//If an event is created in a cold region or outside the physical volume, then it will be killed.
	if((volumeName != "phantom_physicalV" || volumeName == "coldRegion_physicalV") && particleName== "e+"){
		track->SetTrackStatus(fStopAndKill);
		return;		
	}

	//Get the source position if the process name is annihilation
	if(processName == "annihil" && volumeName=="phantom_physicalV"){
		doiPETAnalysis::GetInstance()->SetSourcePosition(pos);
		run->SetAnnihilationTime((track->GetGlobalTime()));
	}

	//get event ID
	//doiPETAnalysis::GetInstance()->SetEventID (eventID);

	//Get scatter information in the phantom by the annihilation photon before detected by the detector. Note that the scatter index is initialized to 0. 
	//If there is scatter, the index is 1, and if not it is 0.
	if(edep>0 && (volumeName == "phantom_physicalV")) scatterIndex = 1;
	doiPETAnalysis::GetInstance()->GetScatterIndexInPhantom(scatterIndex);

	///////////////////    Retrive (Extract) information in the crystal ///////////////////////////
	if(edep>0. && volumeName=="Crystal_physicalV"){

		//get the copy number of the block
		blockID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(2);

		//get the crystal copy number
		crystalID = aStep->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber(0);

		//G4cout<<localPosition.x()<<" "<<localPosition.y()<<" "<<localPosition.z()<<G4endl;

		//get event ID
		//doiPETAnalysis::GetInstance()->GetEventID (eventID);

		//define a pointer for extracting information 
		InteractionInformation* ExtractIntInfo = new InteractionInformation();

		//get the energy deposition of each interaction
		ExtractIntInfo->SetEdep( edep );


		//get the blockID
		ExtractIntInfo->SetBlockNo( blockID );


		//get the crystal ID
		ExtractIntInfo->SetCrystalNo( crystalID );


		//get local position interaction in the crystal
		ExtractIntInfo->SetInteractionPositionInCrystal(localPosition);

		//get the global time of the interaction with the crystal
		ExtractIntInfo->SetGlobalTime( track->GetGlobalTime()); 
		//G4cout<<"Step: "<<eventID<<" "<<blockID<<" "<<crystalID<<" "<<edep<<G4endl;

		//pass all the obtained information to the doiPETAnalysis class
		//doiPETAnalysis::GetInstance()->GetIntractionInfomation(ExtractIntInfo);
		run->GetIntractionInfomation(ExtractIntInfo);

	}
}
