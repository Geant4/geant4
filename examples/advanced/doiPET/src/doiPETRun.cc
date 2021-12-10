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
/// \file electromagnetic/TestEm11/src/doiPETRun.cc
/// \brief Implementation of the doiPETRun class
//
// $Id: doiPETRun.cc 71376 2013-06-14 07:44:50Z maire $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "doiPETRun.hh"
#include "doiPETDetectorConstruction.hh"
#include "doiPETPrimaryGeneratorAction.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "doiPETAnalysis.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

doiPETRun::doiPETRun()
	: G4Run()
{
	//
	totalTime = 0 * s;
	prev_totalTime = 0 * s;
	totalEdep = 0*keV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

doiPETRun::~doiPETRun()
{}

//
void doiPETRun::GetIntractionInfomation(InteractionInformation*step)
{
	//The photon can have many interactions in one block. This stores all the interaction (including multiple copies of the same block.)
	//Save the interacting block in the associative container including multiple copies. Since there will be multiple interactions in the crystals wthin a single block
	mapBlockInteraction.insert( std::multimap< G4int, InteractionInformation* >::value_type( step->GetBlockNo(),  step) );

	//The following stores the block number uniquely. Save the interacting block in the associative container without saving multiple copies
	setBlockInteraction.insert( step->GetBlockNo());
}

///////// InteractingCrystal ///////////////////////////////////

void doiPETRun::FindInteractingCrystal()
{
	std::multimap< G4int, InteractionInformation* >::iterator mitr;
	std::set<G4int>::iterator sitr;
	G4int cystalIDTemp;

	for(sitr=setBlockInteraction.begin(); sitr!=setBlockInteraction.end(); sitr++){
		blockID = (*sitr);
		mitr = mapBlockInteraction.find(blockID);
		numberofInteractions = mapBlockInteraction.count(blockID); //number of interactions in which the interaction deposite energy in a block.
		totalEdep=edepMax=0.;
		for(G4int i=0; i<numberofInteractions; i++, mitr++){
			edep = (*mitr).second->GetEdep();
			cystalIDTemp = (*mitr).second->GetCrystalNo();
			edep_AfterCrystalBlurring = doiPETAnalysis::GetInstance()->QuantumEffifciency(edep, blockID, cystalIDTemp);
			totalEdep += edep_AfterCrystalBlurring;

			crystalID_vec.push_back(cystalIDTemp);
			posInter_vec.push_back(((*mitr).second->GetInteractionPositionInCrystal()));
			edepInCry_vec.push_back(edep_AfterCrystalBlurring);

			//This is  to identify the crystal with the highest energy deposition. 
			if(edepMax<edep_AfterCrystalBlurring){
				edepMax=edep_AfterCrystalBlurring;
				crystalID = cystalIDTemp;
				interactionTime = (*mitr).second->GetGlobalTime();

				//position of interaction based on highest edep. See below for energy weighted position of interaction
				//interactionPos = (*mitr).second->GetInteractionPositionInCrystal();
			}
		}
		if(totalEdep==0)continue;
		
		//The interaction position is based on energy weighted center of mass calculation
		interactionPos = CenterOfMassInteractionPos(crystalID_vec, edepInCry_vec, totalEdep, posInter_vec);
		doiPETAnalysis::GetInstance()->ReadOut(blockID,crystalID,interactionTime,timeOfAnnihil,interactionPos,totalEdep);

		crystalID_vec.clear();
		posInter_vec.clear();
		edepInCry_vec.clear();
	}
	
}

//Apply center of mass calculation in determining the interaction position
G4ThreeVector doiPETRun::CenterOfMassInteractionPos(const std::vector<G4int> &cryID, const std::vector<G4double> &energyDep, G4double edepTotal, const std::vector<G4ThreeVector> &pos){
	posInterInCrystal = G4ThreeVector();
	for(std::size_t i=0; i<cryID.size();i++){
		posInterInCrystal += pos[i] * energyDep[i];
	}
	return (posInterInCrystal)/edepTotal;
}


void doiPETRun::SetAnnihilationTime(G4double t){
	timeOfAnnihil = t;
}
//
void doiPETRun::Clear(){
	//clear the content of the associative containers, the multimap and the set
	std::multimap< G4int, InteractionInformation* >::iterator mitr2;
	for( mitr2= mapBlockInteraction.begin(); mitr2!=mapBlockInteraction.end(); mitr2++){
		delete (*mitr2).second;
	}

	mapBlockInteraction.clear();
	setBlockInteraction.clear();
	//
}
//
void doiPETRun::GetEventIDRun(G4int evt){
	eventID = evt;
}
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void doiPETRun::Merge(const G4Run* run)
{	
	G4Run::Merge(run); 
} 
