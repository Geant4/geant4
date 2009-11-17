#ifndef CML2DummySDH
#define CML2DummySDH

#include "G4VSensitiveDetector.hh"
#include "G4TouchableHistory.hh"
#include "G4VTouchable.hh"

class CML2DummySD : public G4VSensitiveDetector
{
public:
	CML2DummySD(G4String name="");
	~CML2DummySD(void);
	G4bool ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist);
};

#endif
