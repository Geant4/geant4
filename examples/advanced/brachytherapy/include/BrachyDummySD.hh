//    ********************************
//    *                              *
//    *       BrachyDummySD.hh       *
//    *                              *
//    ********************************

// Dummy sensitive used only to flag sensitivity in cells of RO geometry.

#ifndef BrachyDummySD_h
#define BrachyDummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class BrachyDummySD : public G4VSensitiveDetector
{
 public:
 	BrachyDummySD();
  	~BrachyDummySD() {};
  
	void Initialize(G4HCofThisEvent*HCE) {};
	G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
	void EndOfEvent(G4HCofThisEvent*HCE) {};
	void clear() {};
	void DrawAll() {};
	void PrintAll() {};
};

BrachyDummySD::BrachyDummySD() : G4VSensitiveDetector("dummySD")
{}
#endif
