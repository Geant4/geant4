#ifndef DetectorSliceCalorimeterHit_h
#define DetectorSliceCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"

class DetectorSliceCalorimeterHit; 
typedef G4THitsCollection<DetectorSliceCalorimeterHit> DetectorSliceCalorimeterHitsCollection;


class DetectorSliceCalorimeterHit : public G4VHit {

public:

  DetectorSliceCalorimeterHit();
  ~DetectorSliceCalorimeterHit();
  DetectorSliceCalorimeterHit( const DetectorSliceCalorimeterHit &right );
  // Constructors and Destructor.

  const DetectorSliceCalorimeterHit& operator=( const DetectorSliceCalorimeterHit & right );
  // Assignment operator.

public:
  
  void Draw();
  void Print();
  
  inline void SetEdep( const G4double de );
  inline void AddEdep( const G4double de );
  inline G4double GetEdep() const;
  // Set/Add/Get methods for the deposited energy.

private:

  G4double edep;  // Sum of all deposited energies.  

};


inline void DetectorSliceCalorimeterHit::SetEdep(const G4double de) { 
  edep = de;
}

inline void DetectorSliceCalorimeterHit::AddEdep(const G4double de) { 
  edep += de;
}

inline G4double DetectorSliceCalorimeterHit::GetEdep() const { 
  return edep;
}


#endif


