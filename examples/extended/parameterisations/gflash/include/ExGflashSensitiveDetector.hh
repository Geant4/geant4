#ifndef TILE05SENSITIVEDETECTOR_H
#define TILE05SENSITIVEDETECTOR_H

using namespace std;
#include "G4VSensitiveDetector.hh"
#include "G4VGFlashSensitiveDetector.hh"
#include "G4GFlashSpot.hh"
#include "ExGflashDetectorConstruction.hh"
#include "ExGflashHitsCollection.hh"
#include "globals.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class ExGflashSensitiveDetector: public G4VSensitiveDetector ,public G4VGFlashSensitiveDetector{
public:
	ExGflashSensitiveDetector(G4String, ExGflashDetectorConstruction* det);
	~ExGflashSensitiveDetector();
	
	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	G4bool ProcessHits(G4GFlashSpot*aSpot,G4TouchableHistory*); 
	void EndOfEvent(G4HCofThisEvent*);
private:
	ExGflashHitsCollection* caloHitsCollection;
	G4double edep;
	ExGflashDetectorConstruction* Detector;
};
#endif
