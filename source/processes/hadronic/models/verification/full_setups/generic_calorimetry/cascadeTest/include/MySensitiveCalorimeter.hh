#ifndef MySensitiveCalorimeter_h
#define MySensitiveCalorimeter_h 1

#include "G4VSensitiveDetector.hh"
#include "MyCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class MySensitiveCalorimeter : public G4VSensitiveDetector {

public:

  MySensitiveCalorimeter(G4String name);
  ~MySensitiveCalorimeter();

  void Initialize(G4HCofThisEvent* HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
  void EndOfEvent(G4HCofThisEvent* HCE);

  void clear();
  void DrawAll();
  void PrintAll();

private:

  MyCalorimeterHitsCollection* calCollection;

  bool hitExists;

};


#endif
