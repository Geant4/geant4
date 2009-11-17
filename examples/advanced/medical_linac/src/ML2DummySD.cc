#include "ML2DummySD.h"

CML2DummySD::CML2DummySD(G4String name) : G4VSensitiveDetector(name)
{
}

CML2DummySD::~CML2DummySD(void)
{
}
G4bool CML2DummySD::ProcessHits(G4Step *aStep, G4TouchableHistory *ROHist)
{return true;}
