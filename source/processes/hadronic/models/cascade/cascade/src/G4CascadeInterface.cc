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
						       G4V3DNucleus* theNucleus) 
{
  return NULL;
}

G4VParticleChange* G4CascadeInterface::
ApplyYourself(const G4Track& aTrack, G4Nucleus& theNucleus)
{
  if (verboseLevel > 3) 
  {
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;
  }

  // We make conversion between native Geant4 and Bertini cascade classes.
  theResult.Initialize(aTrack);

  enum particleType { proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7 };
  G4int bulletType = 0;

    // Coding particles 
    if(aTrack.GetDefinition() ==    G4Proton::Proton()    ) bulletType = proton;
    if(aTrack.GetDefinition() ==   G4Neutron::Neutron()   ) bulletType = neutron;
    if(aTrack.GetDefinition() ==  G4PionPlus::PionPlus()  ) bulletType = pionPlus;
    if(aTrack.GetDefinition() == G4PionMinus::PionMinus() ) bulletType = pionMinus;
    if(aTrack.GetDefinition() ==  G4PionZero::PionZero()  ) bulletType = pionZero;

    // Code momentum and energy 
    // NOTE: Geant4 units are MeV = 1 and GeV = 1000. Cascade codes by default use MeV.
    G4std::vector<G4double> momentumBullet(4);
    momentumBullet[0] =aTrack.GetDynamicParticle()->Get4Momentum().e();
    momentumBullet[1] =aTrack.GetDynamicParticle()->Get4Momentum().px();
    momentumBullet[2] =aTrack.GetDynamicParticle()->Get4Momentum().py();
    momentumBullet[3] =aTrack.GetDynamicParticle()->Get4Momentum().pz();

    G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, bulletType); 

    // Set target
    G4std::vector<G4double> targetMomentum(4, 0.0);

    G4InuclNuclei * target = new G4InuclNuclei(targetMomentum, 
					       theNucleus.GetN(), 
        				       theNucleus.GetZ());
    target->setEnergy();

    // Resigister collider
    G4ElementaryParticleCollider* collider = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*        cascader = new G4IntraNucleiCascader;
 
    cascader->setElementaryParticleCollider(collider);
    cascader->setInteractionCase(1); // Interaction type is particle with nuclei.

    // Make INC
    G4CollisionOutput output =  cascader->collide(bullet, target); 

    if (verboseLevel > 1) 
    {
      G4cout << " Bertini cascade: " << G4endl;
      output.printCollisionOutput();
    }

    // Convert Bertini data to Geant4 format
    G4std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();

    if(!particles.empty()) 
    { 
      particleIterator ipart;
      G4int outgoingParticle;

      theResult.SetNumberOfSecondaries(particles.size());
      for(ipart = particles.begin(); ipart != particles.end(); ipart++) 
      {
	outgoingParticle = ipart->type();
	G4std::vector<G4double> mom = ipart->getMomentum();
	G4double ekin = ipart->getKineticEnergy();

	G4DynamicParticle* cascadeParticle = NULL;

	switch(outgoingParticle)
	  {
	  case proton: 
	    cascadeParticle = 
	      new G4DynamicParticle(G4Proton::ProtonDefinition(), 
				    G4ParticleMomentum(mom[1], mom[2], mom[3]), ekin);
	    break; 
	  case neutron: 
	    cascadeParticle = 
	      new G4DynamicParticle(G4Neutron::NeutronDefinition(), 
				    G4ParticleMomentum(mom[1], mom[2], mom[3]), ekin);
	    break;
	  case pionPlus: 
	    cascadeParticle = 
	      new G4DynamicParticle(G4PionPlus::PionPlusDefinition(), 
				    G4ParticleMomentum(mom[1], mom[2], mom[3]), ekin);
	    break;
	  case pionMinus:
	    cascadeParticle = 
	      new G4DynamicParticle(G4PionMinus::PionMinusDefinition(), 
				    G4ParticleMomentum(mom[1], mom[2], mom[3]), ekin);
	    break;
	  case pionZero: 
	    cascadeParticle = 
	      new G4DynamicParticle(G4PionZero::PionZeroDefinition(), 
				    G4ParticleMomentum(mom[1], mom[2], mom[3]), ekin);
	    break;
	  default: cout << " ERROR: G4CascadeInterface::Propagate undefined particle type";
	  }

	theResult.AddSecondary(cascadeParticle); 
      }
    }

  return &theResult;
}

