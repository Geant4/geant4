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

typedef vector<G4InuclElementaryParticle>::iterator particleIterator;
      
G4VParticleChange* G4CascadeInterface::ApplyYourself(
						     const G4Track& aTrack, G4Nucleus& theNucleus) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;
  }

  G4std::cout << "Please remove from your physics list."<<G4endl;
  G4Exception("SEVERE: G4CascadeInterface model interface called stand-allone.");
  return new G4ParticleChange;
}
   
G4ReactionProductVector* G4CascadeInterface::Propagate(G4KineticTrackVector* theSecondaries, G4V3DNucleus* theNucleus) {



 if (verboseLevel > 3) {
    G4cout << " >>> G4CascadeInterface::Propagate" << G4endl;
  }

  G4ReactionProductVector * theTotalResult = new G4ReactionProductVector;
 
  // Set bullet
  //      1 - proton
  //      2 - neutron
  //      3 - pi+
  //      5 - pi-
  //      7 - pi 0
  G4int particleType;

  for (G4int list = 0; list < theSecondaries->entries(); list++) {
    G4KineticTrack *aTrack = theSecondaries->at(list);
    // Type 
    if(aTrack->GetDefinition() == G4Proton::Proton()) particleType=1;
    if(aTrack->GetDefinition() == G4Neutron::Neutron()) particleType=2;
    G4InuclElementaryParticle particle(particleType);

    // Momentum 
  
    //  theNew->SetTotalEnergy(aTrack->Get4Momentum().e());
    vector<G4double> momentumBullet(4);
    momentumBullet[0] =aTrack->Get4Momentum().e();
    momentumBullet[1] =aTrack->Get4Momentum().px();
    momentumBullet[2] =aTrack->Get4Momentum().py();
    momentumBullet[3] =aTrack->Get4Momentum().pz();

    G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, particleType); 

    // Interface to Bertini cascade converts between G4 and cascade classes.

    // Set target

    vector<G4double> targetMomentum(4, 0.0);
    //    targetMomentum[0] = theNucleus->Get4Momentum().e();

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

    vector<G4InuclElementaryParticle> particles = 
      output.getOutgoingParticles();

    // Convert Bertini data to G4 internal format
    if(!particles.empty()) { 
      particleIterator ipart;
      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {

	vector<G4double> mom = ipart->getMomentum();
	G4double ekin = ipart->getKineticEnergy();

	// if particle type is neutron
	G4DynamicParticle * aNeutron = new G4DynamicParticle(G4Neutron::NeutronDefinition(),
					 G4ParticleMomentum(mom[1], mom[2], mom[3]),
					 ekin*MeV);
	theParticleChange.AddSecondary(aNeutron);
      };
    };
  };
  // Fill cascade part into the result, and return
  //  for(G4int ll = 0; ll < aPreResult->entries(); ll++) {
  //    theTotalResult->insert(aPreResult->at(ll));
//  }

  return theTotalResult;
}



/*

  // Decay the strong resonances
  G4KineticTrackVector *result1, *secondaries, *result;
  result1 = theSecondaries;
  result = new G4KineticTrackVector();
     
  for (G4int aResult = 0; aResult < result1->entries(); aResult++) {
  G4ParticleDefinition * pdef;
  pdef=result1->at(aResult)->GetDefinition();
  secondaries = NULL;
  if ( pdef->IsShortLived() ) {
  secondaries = result1->at(aResult)->Decay(); // Make decay
  }
  if ( secondaries == NULL ) {
  result->insert(result1->at(aResult));
  result1->at(aResult) = NULL; // Protect for clearAndDestroy 
  } else {
  for (G4int aSecondary = 0; aSecondary < secondaries->entries(); aSecondary++) {
  result1->append(secondaries->at(aSecondary));
  }
  delete secondaries;
  }
  }
  result1->clearAndDestroy();
  delete result1;
     
  // Prepare the fragment
  G4Fragment anInitialState;
  G4int anA = theNucleus->GetMassNumber();
  G4int aZ = theNucleus->GetCharge();
  G4int numberOfEx = 0;
  G4int numberOfCh = 0;
  G4int numberOfHoles = 0;
  G4double exEnergy = 0;
  G4ThreeVector exciton3Momentum(0, 0, 0);

  // Loop over secondaries
  for(G4int list = 0; list < result->entries(); list++) {
  G4KineticTrack *aTrack = result->at(list);
  if(aTrack->GetDefinition() != G4Proton::Proton() && 
  aTrack->GetDefinition() != G4Neutron::Neutron()) {
  G4ReactionProduct * theNew = new G4ReactionProduct(aTrack->GetDefinition());
  theNew->SetMomentum(aTrack->Get4Momentum().vect());
  theNew->SetTotalEnergy(aTrack->Get4Momentum().e());
  theTotalResult->insert(theNew);            
  } else if(aTrack->Get4Momentum().t() - aTrack->Get4Momentum().mag() > 80 * MeV) {
  G4ReactionProduct * theNew = new G4ReactionProduct(aTrack->GetDefinition());
  theNew->SetMomentum(aTrack->Get4Momentum().vect());
  theNew->SetTotalEnergy(aTrack->Get4Momentum().e());
  theTotalResult->insert(theNew);            
  } else if(aTrack->GetPosition().mag() > theNucleus->GetNuclearRadius()) {
  G4ReactionProduct * theNew = new G4ReactionProduct(aTrack->GetDefinition());
  theNew->SetMomentum(aTrack->Get4Momentum().vect());
  theNew->SetTotalEnergy(aTrack->Get4Momentum().e());
  theTotalResult->insert(theNew);            
  } else {
  // Within the nucleus, neutron or proton
  // now calculate  A, Z of the fragment, momentum, number of exciton states
  anA++;;
  numberOfEx++;
  aZ += aTrack->GetDefinition()->GetPDGCharge();
  numberOfCh += aTrack->GetDefinition()->GetPDGCharge();
  exciton3Momentum += aTrack->Get4Momentum().vect();
  exEnergy += (aTrack->Get4Momentum().t()-aTrack->Get4Momentum().m());
  }
  }
     
  // Loop over wounded nucleus
  G4Nucleon * theCurrentNucleon = theNucleus->StartLoop() ? theNucleus->GetNextNucleon() : NULL;
  while(theCurrentNucleon != NULL) {
  if(theCurrentNucleon->AreYouHit()) {
  numberOfHoles++;
  numberOfEx++;
  anA--;
  aZ -= theCurrentNucleon->GetDefinition()->GetPDGCharge();
  exciton3Momentum -= theCurrentNucleon->Get4Momentum().vect();
  exEnergy+=theCurrentNucleon->GetBindingEnergy();
  }
  theCurrentNucleon = theNucleus->GetNextNucleon();
  }   
     
  G4double residualMass =  
  G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(aZ, anA);
  residualMass += exEnergy;
  G4LorentzVector exciton4Momenum(exciton3Momentum, 
  sqrt(exciton3Momentum.mag2() + residualMass * residualMass));
     
  anInitialState.SetA(anA);
  anInitialState.SetZ(aZ);
  anInitialState.SetNumberOfCharged(numberOfCh);
  anInitialState.SetNumberOfHoles(numberOfHoles);
  anInitialState.SetNumberOfExcitons(numberOfEx);
  anInitialState.SetMomentum(exciton4Momentum);
  // anInitialState.SetExcitationEnergy(exEnergy); // now a redundant call.

  // Call pre-compound
  const G4Fragment aFragment(anInitialState);
  if(theDeExcitation) {
  G4ReactionProductVector * aPreResult = theDeExcitation->DeExcite(aFragment);
  // Fill pre-compound part into the result, and return
  for(G4int ll = 0; ll < aPreResult->entries(); ll++) {
  theTotalResult->insert(aPreResult->at(ll));
  }
  delete aPreResult;
  } else {
  // G4Exception("Please register an evaporation phase with G4CascadeInterface.");
  }
  // now return
    
  result->clearAndDestroy();
  delete result;

*/
