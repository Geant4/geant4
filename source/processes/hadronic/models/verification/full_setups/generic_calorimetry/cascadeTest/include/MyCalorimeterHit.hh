#ifndef MyCalorimeterHit_h
#define MyCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"

class MyCalorimeterHit; 
typedef G4THitsCollection<MyCalorimeterHit> MyCalorimeterHitsCollection;


class MyCalorimeterHit : public G4VHit {

public:

  MyCalorimeterHit();
  ~MyCalorimeterHit();
  MyCalorimeterHit(const MyCalorimeterHit &right);
  // Constructors and Destructor.

  const MyCalorimeterHit& operator=(const MyCalorimeterHit &right);
  // Assignment operator.

public:
  
  void Draw();
  void Print();
  
  inline void SetEdep(const G4double de);
  inline void AddEdep(const G4double de);
  inline G4double GetEdep() const;
  // Set/Add/Get methods for the deposited energy.

private:

  G4double edep;  // Sum of all deposited energies.  

};


inline void MyCalorimeterHit::SetEdep(const G4double de) { 
  edep = de;
}

inline void MyCalorimeterHit::AddEdep(const G4double de) { 
  edep += de;
}

inline G4double MyCalorimeterHit::GetEdep() const { 
  return edep;
}

#endif


