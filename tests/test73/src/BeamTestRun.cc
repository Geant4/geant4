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
	G4cout << "Wow you did it :) ... Have a beer!" << G4endl;


}


