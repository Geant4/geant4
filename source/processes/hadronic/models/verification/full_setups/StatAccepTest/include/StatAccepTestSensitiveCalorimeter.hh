#ifndef StatAccepTestSensitiveCalorimeter_h
#define StatAccepTestSensitiveCalorimeter_h 1

#include "G4VSensitiveDetector.hh"
#include "StatAccepTestCalorimeterHit.hh"

class G4Step;
class G4HCofThisEvent;
class G4TouchableHistory;


class StatAccepTestSensitiveCalorimeter : public G4VSensitiveDetector {

public:

  StatAccepTestSensitiveCalorimeter( const G4String name, 
				     const G4int numberOfReplicasIn );
  ~StatAccepTestSensitiveCalorimeter();

  void Initialize( G4HCofThisEvent* HCE );
  G4bool ProcessHits( G4Step* aStep, G4TouchableHistory* ROhist );
  void EndOfEvent( G4HCofThisEvent* HCE );

  void clear();
  void DrawAll();
  void PrintAll();

  inline void Update( const G4int numberOfReplicasIn );
  // Update the change of the geometry.

private:

  StatAccepTestCalorimeterHitsCollection* calCollection;
  G4int numberOfReplicas;
  std::vector< G4int > positionID;

};


inline void StatAccepTestSensitiveCalorimeter::Update( const G4int numberOfReplicasIn ) {
  numberOfReplicas = numberOfReplicasIn;
}

#endif
