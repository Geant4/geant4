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
/// \file electromagnetic/TestEm11/include/doiPETRun.hh
/// \brief Definition of the doiPETRun class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef doiPETRun_h
#define doiPETRun_h 1

#include "G4Run.hh"
#include "G4VProcess.hh"
#include "globals.hh"

#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <set>
class doiPETDetectorConstruction;
class G4ParticleDefinition;
class doiPETRunAction;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class InteractionInformation; 
class doiPETAnalysis;
class doiPETRun : public G4Run
{
public:
	doiPETRun(/*doiPETDetectorConstruction**/);
	~doiPETRun();

public:
	//
	void GetIntractionInfomation(InteractionInformation*);
	void FindInteractingCrystal();
	void Clear();
	void OpenRun(G4String);
	void GetEventIDRun(G4int);
	void CalulateAcquisitionTime();
	void SetAnnihilationTime(G4double);
	void SetEventID(G4int);
	virtual void Merge(const G4Run*);
	// void EndOfRun();     

private:

	std::multimap< G4int, InteractionInformation* > mapBlockInteraction;
	std::set<G4int> setBlockInteraction;
	std::ofstream ofs;
	G4int eventID;
	doiPETRunAction* fRunAction;
	doiPETAnalysis* fAnalysis;
	G4double totalEdep;
	G4int blockID, crystalID;

	//
	G4double activityNow;
	G4double InitialActivity;
	G4double halfLife;

	G4double totalTime;
	G4double prev_totalTime;
	G4double timeInterval;

	G4ThreeVector interactionPos;
	//
	G4double interactionTime;
	G4double timeOfAnnihil;

	G4int numberofInteractions;
	G4double edep;
	G4double edepMax;
	G4double edep_AfterCrystalBlurring;

};


//The following is to get interaction information
class InteractionInformation
{

public:
	InteractionInformation(){;};
	~InteractionInformation(){;};

	//set energy deposition
	void SetEdep(G4double e) {edep=e;};

	//set block number
	void SetBlockNo(G4int n) {blockNo=n;};

	//set crystal ID (this is continuous numbering of crystals)
	void SetCrystalNo(G4int n) {crystalNo=n;};

	//set global time
	void SetGlobalTime(G4double t){globalTime=t;};

	//set interaction position with respect to the crystal axis (local position)
	void SetInteractionPositionInCrystal(G4ThreeVector pos){crystalPosition = pos;};


	G4double GetEdep() {return edep;};
	G4int GetBlockNo() {return blockNo;};
	G4int GetCrystalNo() {return crystalNo;};
	G4double GetGlobalTime() {return globalTime;};

	//Interaction position in the detector
	G4ThreeVector GetInteractionPositionInCrystal(){return crystalPosition;};


private:
	G4double edep;
	G4double globalTime;
	G4int blockNo;
	G4int crystalNo;
	G4ThreeVector crystalPosition;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

