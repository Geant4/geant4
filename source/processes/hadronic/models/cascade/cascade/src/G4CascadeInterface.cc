#include "G4CascadeInterface.hh"

#include "globals.hh"
#include "G4DynamicParticleVector.hh"
#include "G4IonTable.hh"

#include "G4IntraNucleiCascader.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4V3DNucleus.hh"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4CascadeInterface::G4CascadeInterface()
  :verboseLevel(1)  {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::G4CascadeInterface" << G4endl;
  }
};
   
G4ReactionProductVector* G4CascadeInterface::Propagate(G4KineticTrackVector* theSecondaries, 
						       G4V3DNucleus* theNucleus) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::Propagate" << G4endl;
  }

  // We make conversion between native Geant4 and Bertini cascade classes.

  enum particleType { proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7 };

  G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;

  G4int particleType;

  for (G4int list = 0; list < theSecondaries->entries(); list++) {
    G4KineticTrack *aTrack = theSecondaries->at(list);

    // Coding particles 
    if(aTrack->GetDefinition() ==    G4Proton::Proton())    particleType = proton;
    if(aTrack->GetDefinition() ==   G4Neutron::Neutron())   particleType = neutron;
    if(aTrack->GetDefinition() ==  G4PionPlus::PionPlus())  particleType = pionPlus;
    if(aTrack->GetDefinition() == G4PionMinus::PionMinus()) particleType = pionMinus;
    if(aTrack->GetDefinition() ==  G4PionZero::PionZero())  particleType = pionZero;

    G4InuclElementaryParticle particle(particleType);

    // Code momentum   
    G4std::vector<G4double> momentumBullet(4);
    momentumBullet[0] =aTrack->Get4Momentum().e();
    momentumBullet[1] =aTrack->Get4Momentum().px();
    momentumBullet[2] =aTrack->Get4Momentum().py();
    momentumBullet[3] =aTrack->Get4Momentum().pz();

    G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, particleType); 

    // Set target
    G4std::vector<G4double> targetMomentum(4, 0.0);
    //targetMomentum[0] = theNucleus->Get4Momentum().e();

    G4InuclNuclei * target = new G4InuclNuclei(
					       targetMomentum, 
					       theNucleus->GetMassNumber(), 
					       theNucleus->GetCharge());
    target->setEnergy();
    G4ElementaryParticleCollider* collider = new G4ElementaryParticleCollider;

    // Resigister collider
    G4IntraNucleiCascader *  cascader = new G4IntraNucleiCascader; 
    cascader->setElementaryParticleCollider(collider);
    cascader->setInteractionCase(1); // Interaction type is particle with nuclei.

    // Actual cascade generation 
    G4CollisionOutput output =  cascader->collide(bullet, target); 

    if (verboseLevel > 1) {
      G4cout << " Bertini cascade: " << G4endl;
      output.printCollisionOutput();
    }

    G4std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();

    // Convert Bertini data to Geant4 format
    if(!particles.empty()) { 
      particleIterator ipart;
      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {

	G4std::vector<G4double> mom = ipart->getMomentum();
	G4double ekin = ipart->getKineticEnergy();

	// Particle type is neutron
	G4DynamicParticle * aNeutron = new G4DynamicParticle(G4Neutron::NeutronDefinition(),
							     G4ParticleMomentum(mom[1], mom[2], mom[3]),
							     ekin * MeV);
	theParticleChange.AddSecondary(aNeutron);
      };
    };
  };

  // Fill cascade part into the result, and return
  //  for(G4int ll = 0; ll < aPreResult->entries(); ll++) {
  //    theTotalResult->insert(aPreResult->at(ll)); }

  return theTotalResult;
}

G4VParticleChange* G4CascadeInterface::ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;
  }

  G4std::cout << "Please remove from your physics list."<<G4endl;
  G4Exception("SEVERE: G4CascadeInterface model interface called stand-allone.");

  return new G4ParticleChange;
}
