//
// Dummy sensitive used only to flag sensitivity
// in cells of RO geometry.
//

#ifndef ExN04DummySD_h
#define ExN04DummySD_h 1

#include "G4VSensitiveDetector.hh"
class G4Step;

class ExN04DummySD : public G4VSensitiveDetector
{
public:
  ExN04DummySD();
  ~ExN04DummySD() {};
  
  void Initialize(G4HCofThisEvent*HCE) {};
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
  void EndOfEvent(G4HCofThisEvent*HCE) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};
};
ExN04DummySD::ExN04DummySD()
  : G4VSensitiveDetector("dummySD")
{}
#endif
