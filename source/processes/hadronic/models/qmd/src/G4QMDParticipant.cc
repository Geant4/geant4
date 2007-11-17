#include "G4QMDParticipant.hh"

G4QMDParticipant::G4QMDParticipant( G4ParticleDefinition* pd , G4ThreeVector p , G4ThreeVector r )
: definition ( pd )
, momentum ( p )
, position ( r )
, projectile ( false )
, target ( false )
{
   ; 
}



G4QMDParticipant::~G4QMDParticipant()
{
   ;
}



G4LorentzVector G4QMDParticipant::Get4Momentum()
{
   G4LorentzVector p4 ( momentum , std::sqrt ( std::pow ( definition->GetPDGMass()/GeV , 2 ) + momentum*momentum ) );
   return p4;
}
