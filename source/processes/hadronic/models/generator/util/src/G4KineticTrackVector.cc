#include "G4KineticTrackVector.hh"

G4KineticTrackVector::G4KineticTrackVector()
    {
    }


//****************************************************************************************************************
// These methods were implemented by Maxim Komogorov
// Maxim.Komogorov@cern.ch

void G4KineticTrackVector::BoostBeam(G4ThreeVector& BeamMom)
    {
    for(G4int c1 = 0; c1 < length(); c1++)
        {
        G4KineticTrack& KT =*at(c1);
        G4LorentzVector Mom = KT.Get4Momentum();        
        G4ThreeVector Velocity = (1/sqrt(BeamMom.mag2() + sqr(KT.GetDefinition()->GetPDGMass())))*BeamMom;
        Mom.boost(Velocity);
        KT.Set4Momentum(Mom);
        }
    }
//--------------------------------------------------------------------------------------------------------------

void G4KineticTrackVector::Boost(G4ThreeVector& Velocity)
    { 
    for(G4int c1 = 0; c1 < length(); c1++)
        {
        G4KineticTrack& KT =*at(c1);
        G4LorentzVector Mom = KT.Get4Momentum();        
        Mom.boost(Velocity);
        KT.Set4Momentum(Mom);
        }
    }

//--------------------------------------------------------------------------------------------------------------

void G4KineticTrackVector::Shift(G4ThreeVector& Pos)
    {
    for(G4int c1 = 0; c1 < length(); c1++)
        {
        G4KineticTrack& KT =*at(c1);
        KT.SetPosition(KT.GetPosition() + Pos);
        }
    }
 

//****************************************************************************************************************
