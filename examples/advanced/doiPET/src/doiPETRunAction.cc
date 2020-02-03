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



//
#include "doiPETRunAction.hh"
#include "doiPETRun.hh"
#include "doiPETDetectorConstruction.hh"
#include "doiPETPrimaryGeneratorAction.hh"
#include "doiPETAnalysis.hh"

#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <sstream>
#include <vector>
#include <fstream>

//
doiPETRunAction::doiPETRunAction() : G4UserRunAction()
{}
//


doiPETRunAction::~doiPETRunAction()
{}

//
G4Run* doiPETRunAction::GenerateRun()
{ 
	fRun = new doiPETRun();
	return fRun;
}

//

void doiPETRunAction::BeginOfRunAction(const G4Run* aRun)
{ 
	//inform the runManager to save random number seed
	G4RunManager::GetRunManager()->SetRandomNumberStore(false);
	G4cout << "### Begin of Run " << aRun->GetRunID() << " start." << G4endl;

	if (isMaster) {
		G4Random::showEngineStatus();
	}	
}


void doiPETRunAction::EndOfRunAction(const G4Run*)
{
	G4int nofEvents = fRun->GetNumberOfEvent();
	if (nofEvents == 0) return;
	G4cout<<G4endl<<"### End of Run ("<<nofEvents<<" events)"<<G4endl;
}


