#include "G4BinaryLightIonReaction.hh"
#include "G4LorentzVector.hh"
#include "G4LorentzRotation.hh"
#include <algorithm>


G4VParticleChange *G4BinaryLightIonReaction::
ApplyYourself(const G4Track &aTrack, G4Nucleus & targetNucleus )
{
    // add rotation to Z !!!!! @@@@@@
    
    G4double a1=aTrack.GetDefinition()->GetBaryonNumber();
    G4double z1=aTrack.GetDefinition()->GetPDGCharge();
    G4double m1=aTrack.GetDefinition()->GetPDGMass();
    G4double a2=targetNucleus.GetN();
    G4double z2=targetNucleus.GetZ();
    G4double m2=G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(z2, a2);
    
    G4LorentzVector mom(aTrack.GetDynamicParticle()->Get4Momentum());
    G4LorentzRotation toBreit(mom.boostVector());
        
    G4bool swapped = false;
    if(a2<a1)
    {
      swapped = true;
      G4double tmp(0);
      tmp = a2; a2=a1; a1=tmp;
      tmp = z2; z2=z1; z1=tmp;
      tmp = m2; m2=m1; m1=tmp;
      G4LorentzVector it(m1, G4ThreeVector(0,0,0));
      mom = toBreit*it;
    }

    G4ReactionProductVector * result = 0;
    while(!result)
    {
      G4V3DNucleus * projectile = new G4Fancy3DNucleus;
      projectile->Init(a1, z1);
      G4V3DNucleus * fancyNucleus = new G4Fancy3DNucleus;  
      fancyNucleus->Init(a2, z2);

      G4double impactMax = fancyNucleus->GetOuterRadius()+projectile->GetOuterRadius();
      G4double aX=(2.*G4UniformRand()-1.)*impactMax;
      G4double aY=(2.*G4UniformRand()-1.)*impactMax;
      G4ThreeVector pos(aX, aY, -2.*impactMax);

      G4KineticTrackVector * initalState = new G4KineticTrackVector;
      projectile->DoLorentzBoost(-1.*mom);
      projectile->StartLoop();
      G4Nucleon * aNuc;
      while( (aNuc=projectile->GetNextNucleon()) )
      {
	G4LorentzVector p4 = aNuc->GetMomentum();
	// testing @@@@@
	// p4.setX(0);
	// p4.setY(0);
	// end testing @@@@@
	
	G4KineticTrack * it = new G4KineticTrack(aNuc->GetDefinition(), 0, aNuc->GetPosition()+pos, p4 );
	initalState->push_back(it);
      }
      G4ReactionProductVector *result=theModel.Propagate(initalState, fancyNucleus);
      delete fancyNucleus;
      delete projectile;
      for_each(initalState->begin(), initalState->end(), DeleteKineticTrack());
      delete initalState;
      
      if(result->size()==0) 
      {
        delete result; result=0;
      }
      
    }
    //inverse transformation in case we swapped.
    fancyNucleus->StartLoop();   
    while( (aNuc=fancyNucleus->GetNextNucleon()) )
    {
      if(!aNuc->Hit())
      {
        // stuff them away.
      }
    }
  return 0;
}



