#include "G4PartonPair.hh"

G4PartonPair::G4PartonPair(G4Parton* P1, G4Parton* P2, G4int Type, G4int aDirection)
    {
    CollisionType = Type;
    Parton1 = P1;
    Parton2 = P2;
    Direction = aDirection;
    }

G4PartonPair::G4PartonPair(const G4PartonPair &right)
    {
    G4Exception("You can not make a copy of this object");
    }

G4PartonPair::~G4PartonPair()
    {
    }
 
