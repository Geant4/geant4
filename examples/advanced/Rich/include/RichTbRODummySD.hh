// Rich advanced example for Geant4
// RichTbRODummySD.hh for Rich of LHCb
// History:
// Created: Sajan Easo (Sajan.Easo@cern.ch)
// Revision and changes: Patricia Mendez (Patricia.Mendez@cern.ch)
/////////////////////////////////////////////////////////////////////////////
#ifndef RichTbRODummySD_h
#define RichTbRODummySD_h 1
#include "G4VSensitiveDetector.hh"
class G4Step;

class RichTbRODummySD : public G4VSensitiveDetector
{
public:
  RichTbRODummySD();
  ~RichTbRODummySD() {};

  void Initialize(G4HCofThisEvent*HCE) {};
  G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist) {return false;}
  void EndOfEvent(G4HCofThisEvent*HCE) {};
  void clear() {};
  void DrawAll() {};
  void PrintAll() {};

};
RichTbRODummySD::RichTbRODummySD()
  : G4VSensitiveDetector("RichTbROdummySD") {}

#endif
