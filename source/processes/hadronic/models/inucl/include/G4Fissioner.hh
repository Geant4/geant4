#ifndef G4FISSIONER_HH
#define G4FISSIONER_HH

#include "G4Collider.hh"
#include "G4InuclSpecialFunctions.hh"

using namespace G4InuclSpecialFunctions;

class G4Fissioner : public G4Collider {

public:

G4Fissioner() {};

virtual G4CollisionOutput collide(G4InuclParticle* bullet,
                     G4InuclParticle* target);

private: 

G4double getC2(G4double A1, G4double A2, G4double X3, G4double X4, G4double R12) const; 

G4double getZopt(G4double A1, G4double A2, G4double ZT, 
                   G4double X3, G4double X4, G4double R12) const;
		    
void potentialMinimization(G4double& VP, vector<G4double>& ED, G4double& VC,
   G4double AF, G4double AS, G4double ZF, G4double ZS,
   vector<G4double>& AL1, vector<G4double>& BET1, G4double& R12) const; 

};        

#endif // G4FISSIONER_HH
