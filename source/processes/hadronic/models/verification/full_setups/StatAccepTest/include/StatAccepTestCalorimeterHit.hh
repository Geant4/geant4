#ifndef StatAccepTestCalorimeterHit_h
#define StatAccepTestCalorimeterHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"

class StatAccepTestCalorimeterHit; 
typedef G4THitsCollection<StatAccepTestCalorimeterHit> StatAccepTestCalorimeterHitsCollection;


class StatAccepTestCalorimeterHit : public G4VHit {

public:

  StatAccepTestCalorimeterHit();
  ~StatAccepTestCalorimeterHit();
  StatAccepTestCalorimeterHit( const StatAccepTestCalorimeterHit &right );
  // Constructors and Destructor.

  const StatAccepTestCalorimeterHit& operator=( const StatAccepTestCalorimeterHit & right );
  // Assignment operator.

public:
  
  void Draw();
  void Print();
  
  inline void SetEdep( const G4double de );
  inline void AddEdep( const G4double de );
  inline G4double GetEdep() const;
  // Set/Add/Get methods for the deposited energy.

  inline void SetLayer( const G4int layerNum );
  inline G4int GetLayer() const;
  // Set/Get methods for the layer number. 

private:

  G4double edep;  // Sum of all deposited energies.  
  G4int layer;

};


inline void StatAccepTestCalorimeterHit::SetEdep(const G4double de) { 
  edep = de;
}

inline void StatAccepTestCalorimeterHit::AddEdep(const G4double de) { 
  edep += de;
}

inline G4double StatAccepTestCalorimeterHit::GetEdep() const { 
  return edep;
}


inline void StatAccepTestCalorimeterHit::SetLayer( const G4int layerNum ) { 
  layer = layerNum;
}

inline G4int StatAccepTestCalorimeterHit::GetLayer() const { 
  return layer;
}


#endif


