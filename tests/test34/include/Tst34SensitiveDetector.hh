#ifndef TILE05SENSITIVEDETECTOR_H
#define TILE05SENSITIVEDETECTOR_H

using namespace std;
#include "G4VSensitiveDetector.hh"
#include "G4VGFlashSensitiveDetector.hh"
#include "G4GFlashSpot.hh"
#include "Tst34DetectorConstruction.hh"
#include "Tst34HitsCollection.hh"
#include "globals.hh"
class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;
class Tst34SensitiveDetector: public G4VSensitiveDetector ,public G4VGFlashSensitiveDetector{
public:
	Tst34SensitiveDetector(G4String, Tst34DetectorConstruction* det);
	~Tst34SensitiveDetector();
	
	void Initialize(G4HCofThisEvent*);
	G4bool ProcessHits(G4Step*,G4TouchableHistory*);
	G4bool ProcessHits(G4GFlashSpot*aSpot,G4TouchableHistory*); 
	void EndOfEvent(G4HCofThisEvent*);
private:
	Tst34HitsCollection* caloHitsCollection;
	G4double edep;
	Tst34DetectorConstruction* Detector;
};
#endif
