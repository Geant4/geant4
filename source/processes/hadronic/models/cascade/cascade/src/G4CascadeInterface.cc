#include "G4CascadeInterface.hh"

#include "globals.hh"
#include "G4DynamicParticleVector.hh"
#include "G4IonTable.hh"

#include "G4InuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4V3DNucleus.hh"
#include "G4Track.hh"
#include "G4Nucleus.hh"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4CascadeInterface::G4CascadeInterface()
  :verboseLevel(1)  {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::G4CascadeInterface" << G4endl;
  }
};
   
G4ReactionProductVector* G4CascadeInterface::Propagate(G4KineticTrackVector* theSecondaries, 
						       G4V3DNucleus* theNucleus) {
  return NULL;
};

#define debug_G4CascadeInterface

G4VParticleChange* G4CascadeInterface::ApplyYourself(const G4Track& aTrack, 
						     G4Nucleus& theNucleus) {
#ifdef debug_G4CascadeInterface
  static G4int counter(0);
  counter++;
  G4cerr << "Reaction number "<< counter << G4endl;
#endif
  
  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;
  };

  theResult.Initialize(aTrack);

  // Make conversion between native Geant4 and Bertini cascade classes.
  // NOTE: Geant4 units are MeV = 1 and GeV = 1000. Cascade code by default use GeV = 1.

  enum particleType { proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7 };

  G4int bulletType = 0;

  // Coding particles 
  if(aTrack.GetDefinition() ==    G4Proton::Proton()    ) bulletType = proton;
  if(aTrack.GetDefinition() ==   G4Neutron::Neutron()   ) bulletType = neutron;
  if(aTrack.GetDefinition() ==  G4PionPlus::PionPlus()  ) bulletType = pionPlus;
  if(aTrack.GetDefinition() == G4PionMinus::PionMinus() ) bulletType = pionMinus;
  if(aTrack.GetDefinition() ==  G4PionZero::PionZero()  ) bulletType = pionZero;

  // Code momentum and energy.
  G4std::vector<G4double> momentumBullet(4);
  momentumBullet[0] =aTrack.GetDynamicParticle()->Get4Momentum().e()  / GeV;
  momentumBullet[1] =aTrack.GetDynamicParticle()->Get4Momentum().px() / GeV;
  momentumBullet[2] =aTrack.GetDynamicParticle()->Get4Momentum().py() / GeV;
  momentumBullet[3] =aTrack.GetDynamicParticle()->Get4Momentum().pz() / GeV;

  G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, bulletType); 

  // Set target
  G4std::vector<G4double> targetMomentum(4, 0.0);

  G4InuclNuclei * target = new G4InuclNuclei(targetMomentum, 
					     theNucleus.GetN(), 
					     theNucleus.GetZ());
  target->setEnergy();

  G4CollisionOutput output;

  G4bool runOnlyINC = false; // Flag for running only INC or all models.

  if (runOnlyINC) { // Run only INC
  
    // Resigister collider
    G4ElementaryParticleCollider* colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*     collider = new G4IntraNucleiCascader;
 
    collider->setElementaryParticleCollider(colep);
    collider->setInteractionCase(1); // Interaction type is particle with nuclei.

    output = collider->collide(bullet, target); 
 
  } else { // Run all models

    // Colliders initialisation
    G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*            inc = new G4IntraNucleiCascader; // the actual cascade
    inc->setInteractionCase(1); // Interaction type is particle with nuclei.

    G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
    G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
    G4Fissioner*                     fiss = new G4Fissioner;
    G4BigBanger*                     bigb = new G4BigBanger;

    G4InuclCollider*             collider = new G4InuclCollider(colep, inc, noneq, eqil, fiss, bigb);

    output = collider->collide(bullet, target); 
  }

  if (verboseLevel > 1) {
    G4cout << " Cascade output: " << G4endl;
    output.printCollisionOutput();
  };

  // Convert Bertini data to Geant4 format
  G4std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();

  if(!particles.empty()) { 
    particleIterator ipart;
    G4int outgoingParticle;

    theResult.SetNumberOfSecondaries(particles.size());
    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
      outgoingParticle = ipart->type();
      G4std::vector<G4double> mom = ipart->getMomentum();
      G4double ekin = ipart->getKineticEnergy() * GeV;
      G4ThreeVector aMom(mom[1], mom[2], mom[3]);
      aMom = aMom.unit();

      G4DynamicParticle* cascadeParticle = NULL;

      switch(outgoingParticle) {

      case proton: 
	cascadeParticle = 
	  new G4DynamicParticle(G4Proton::ProtonDefinition(), aMom, ekin);
	break; 

      case neutron: 
	cascadeParticle = 
	  new G4DynamicParticle(G4Neutron::NeutronDefinition(), aMom, ekin);
	break;

      case pionPlus: 
	cascadeParticle = 
	  new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), aMom, ekin);
	break;

      case pionMinus:
	cascadeParticle = 
	  new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), aMom, ekin);
	break;

      case pionZero: 
	cascadeParticle = 
	  new G4DynamicParticle(G4PionZero::PionZeroDefinition(), aMom, ekin);
	break;

      default: cout << " ERROR: G4CascadeInterface::Propagate undefined particle type";
      }

      theResult.AddSecondary(cascadeParticle); 
    }
  }

  return &theResult;
}
