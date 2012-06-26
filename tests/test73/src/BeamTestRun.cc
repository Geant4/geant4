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
// BeamTestRun is a class for accumulating scored quantities over
// an entire run. Event data is accumulated over a run in a G4THitsMap 
// object.
//

#include "BeamTestRun.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4SDManager.hh"
#include "G4VPrimitiveScorer.hh"
#include <assert.h>

BeamTestRun::BeamTestRun(const G4String& detectorName)
{
	// Get the sensitive detector manager
	G4SDManager* manager = G4SDManager::GetSDMpointer();

	// Get the sensitive detector
	G4MultiFunctionalDetector* detector =
		dynamic_cast<G4MultiFunctionalDetector*>(manager->FindSensitiveDetector(detectorName));

	//G4String* detName = dynamic_cast<G4String*>( manager->FindSensitiveDetector(detectorName));
	//G4cout << "Detector Name: " << detName << G4endl;
	G4cout << "Detector Name: " << detector << G4endl;
	// Expect the detector to exist
	//assert (0 != detector);


	// Loop over primitive scorers registered with the detector

		// Get scorer

		// Need to form the full collection name = detector name + "/"+ scorer name 
		// to get the collection id number


		// Expect the collection to have been added to the sensitive detector manager
		// when scorer was registered with G4MultiFunctionalDetector


}

BeamTestRun::~BeamTestRun()
{}


void BeamTestRun::DumpData() const
{
	// Titles
}


