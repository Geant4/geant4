//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
// $Id: G4PreCompoundCascadeInterface.cc,v 1.15 2010-06-17 15:32:35 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100414  M. Kelsey -- Check for K0L/K0S before using G4InuclElemPart::type
// 20100419  M. Kelsey -- Access G4CollisionOutput lists by const-ref, and 
//		const_iterator
// 20100428  M. Kelsey -- Use G4InuclParticleNames enum
// 20100429  M. Kelsey -- Change "case gamma:" to "case photon:"
// 20100517  M. Kelsey -- Follow new ctors for G4*Collider family.
// 20100520  M. Kelsey -- Add missing name string to ctor, follow code changes
//		from G4CascadeInterface.
// 20100617  M. Kelsey -- Rename "debug_" preprocessor flag to G4CASCADE_DEBUG

#include "G4PreCompoundCascadeInterface.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"
#include "G4DynamicParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4LorentzRotation.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4PreCompoundInuclCollider.hh"
#include "G4Track.hh"
#include "G4V3DNucleus.hh"

using namespace G4InuclParticleNames;


typedef std::vector<G4InuclElementaryParticle>::const_iterator particleIterator;
typedef std::vector<G4InuclNuclei>::const_iterator nucleiIterator;

G4PreCompoundCascadeInterface::G4PreCompoundCascadeInterface(const G4String& nam)
  :G4VIntraNuclearTransportModel(nam), verboseLevel(0)  {

  if (verboseLevel > 3) {
    G4cout << " >>> G4PreCompoundCascadeInterface::G4PreCompoundCascadeInterface" << G4endl;
  }
}
   
G4ReactionProductVector* G4PreCompoundCascadeInterface::Propagate(G4KineticTrackVector* , 
								  G4V3DNucleus* ) {
  return 0;
}

G4HadFinalState* 
G4PreCompoundCascadeInterface::ApplyYourself(const G4HadProjectile& aTrack, 
					     G4Nucleus& theNucleus) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4PreCompoundCascadeInterface::ApplyYourself" << G4endl;
  };

#ifdef G4CASCADE_DEBUG_INTERFACE
  static G4int counter(0);
  counter++;
  G4cerr << "Reaction number "<< counter << " "
	 << aTrack.GetDefinition()->GetParticleName() << " "
	 << aTrack.GetKineticEnergy() << G4endl;
#endif

  theResult.Clear();

  G4double eInit     = 0.0;
  G4double eTot      = 0.0;
  G4double sumBaryon = 0.0;
  G4double sumEnergy = 0.0;

  // Make conversion between native Geant4 and Bertini cascade classes.
  // NOTE: Geant4 units are MeV = 1 and GeV = 1000. Cascade code by default use GeV = 1.

  G4int bulletType;
  if (aTrack.GetDefinition() == G4KaonZeroLong::KaonZeroLong() ||
      aTrack.GetDefinition() == G4KaonZeroShort::KaonZeroShort() )
    bulletType = (G4UniformRand() > 0.5) ? kaonZero : kaonZeroBar;
  else 
    bulletType = G4InuclElementaryParticle::type(aTrack.GetDefinition());

  // Code momentum and energy.
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum();
  G4LorentzRotation toZ;
  toZ.rotateZ(-projectileMomentum.phi());
  toZ.rotateY(-projectileMomentum.theta());
  G4LorentzRotation toLabFrame = toZ.inverse();

  G4LorentzVector momentumBullet(0., 0., aTrack.GetTotalMomentum()/GeV,
				 aTrack.GetTotalEnergy()/GeV);

  G4InuclElementaryParticle *  bullet =
    new G4InuclElementaryParticle(momentumBullet, bulletType); 

  sumEnergy = bullet->getKineticEnergy(); // In GeV 
  sumBaryon += bullet->baryon();

  // Set target
  G4InuclNuclei*   target  = 0;
  G4InuclParticle* targetH = 0;

  G4double theNucleusA = theNucleus.GetN();

  if ( !(G4int(theNucleusA) == 1) ) {
    target  = new G4InuclNuclei(theNucleusA, theNucleus.GetZ());
    eInit = bullet->getEnergy() + target->getEnergy();

    sumBaryon += theNucleusA;

    if (verboseLevel > 2) {
      G4cout << "Bullet:  " << G4endl;  
      bullet->printParticle();
    }
    if (verboseLevel > 2) {
      G4cout << "Target:  " << G4endl;  
      target->printParticle();
    }
  }

  G4CollisionOutput output;

  // Colliders initialisation
  G4PreCompoundInuclCollider*  collider = new G4PreCompoundInuclCollider;

  G4int  maxTries = 10; // maximum tries for inelastic collision to avoid infinite loop
  G4int  nTries   = 0;  // try counter

  if (G4int(theNucleusA) == 1) { // special treatment for target H(1,1) (proton)

    targetH = new G4InuclElementaryParticle(1);

    G4float cutElastic[32];

    cutElastic[proton   ]   = 1.0; // 1 GeV
    cutElastic[neutron  ]   = 1.0;
    cutElastic[lambda]      = 1.0;
    cutElastic[sigmaPlus]   = 1.0;
    cutElastic[sigmaZero]   = 1.0;
    cutElastic[sigmaMinus]  = 1.0;
    cutElastic[xiZero]      = 1.0;
    cutElastic[xiMinus]     = 1.0;

    cutElastic[pionPlus ]   = 0.6; // 0.6 GeV

    cutElastic[kaonPlus ]   = 0.5; // 0.5 GeV
    cutElastic[kaonMinus]   = 0.5;
    cutElastic[kaonZero]    = 0.5;
    cutElastic[kaonZeroBar] = 0.5;

    cutElastic[pionMinus]   = 0.2; // 0.2 GeV
    cutElastic[pionZero ]   = 0.2;


    if (momentumBullet.z() > cutElastic[bulletType]) { // inelastic collision possible

      do {   // we try to create inelastic interaction
	output.reset();
	collider->collide(bullet, targetH, output);
	nTries++;
      } while(
	      (nTries < maxTries) &&
	      (output.getOutgoingParticles().size() == 2 && // elastic: bullet + p = H(1,1) coming out
	       (output.getOutgoingParticles().begin()->type() == bulletType || output.getOutgoingParticles().begin()->type() == proton)
	       )
	      );

    } else { // only elastic collision is energetically possible
      collider->collide(bullet, targetH, output);
    }

    sumBaryon += 1;

    eInit = bullet->getEnergy() + target->getEnergy();

    if (verboseLevel > 2) {
      G4cout << "Target:  " << G4endl;
      targetH->printParticle();
    }

  } else {  // treat all other targets excepet H(1,1)

    do  // we try to create inelastic interaction
      {
	output.reset();
	collider->collide(bullet, target, output);
	nTries++;
      } while (
	       (nTries < maxTries) &&
	       output.getOutgoingParticles().size() == 1 &&     // we retry when elastic collision happened
               output.getNucleiFragments().size() == 1 &&            
	       output.getOutgoingParticles().begin()->type() == bullet->type() &&
	       output.getNucleiFragments().begin()->getA() == target->getA() && 
	       output.getNucleiFragments().begin()->getZ() == target->getZ() 
	       );
  }
 
  if (verboseLevel > 1) {
      G4cout << " Cascade output: " << G4endl;
      output.printCollisionOutput();
  }
  
  // Rotate event to put Z axis along original projectile direction
  output.rotateEvent(toLabFrame);

  // Convert cascade data to use hadronics interface
  const std::vector<G4InuclNuclei>& nucleiFragments = output.getNucleiFragments();
  const std::vector<G4InuclElementaryParticle>& particles = output.getOutgoingParticles();

  theResult.SetStatusChange(stopAndKill);

  // Get outcoming particles
  G4DynamicParticle* cascadeParticle = 0;
  if (!particles.empty()) { 
    particleIterator ipart = particles.begin();
    for (; ipart != particles.end(); ipart++) {
      G4int outgoingType = ipart->type();

      eTot += ipart->getEnergy();
      sumBaryon -= ipart->baryon();
      sumEnergy -= ipart->getKineticEnergy();

      if (!ipart->valid() || ipart->quasi_deutron()) {
        G4cerr << " ERROR: G4PreCompoundCascadeInterface::Propagate incompatible"
	       << " particle type " << ipart->type() << G4endl;
	continue;
      }

      // Copy local G4DynPart to public output (handle kaon mixing specially)
      if (outgoingType == kaonZero || outgoingType == kaonZeroBar) {
	G4ThreeVector momDir = ipart->getMomentum().vect().unit();
	G4double ekin = ipart->getKineticEnergy()*GeV;	// Bertini -> G4 units

	G4ParticleDefinition* pd = G4KaonZeroShort::Definition();
	if (G4UniformRand() > 0.5) pd = G4KaonZeroLong::Definition();

	cascadeParticle = new G4DynamicParticle(pd, momDir, ekin);
      } else {
	cascadeParticle = new G4DynamicParticle(ipart->getDynamicParticle());
      }
 
      theResult.AddSecondary(cascadeParticle); 
    }
  }

  // get nuclei fragments
  G4DynamicParticle * aFragment = 0;
  if (!nucleiFragments.empty()) { 
    nucleiIterator ifrag = nucleiFragments.begin();
    for (; ifrag != nucleiFragments.end(); ifrag++) {
      eTot += ifrag->getEnergy();
      sumBaryon -= ifrag->getA();
      sumEnergy -= ifrag->getKineticEnergy();

      if (verboseLevel > 2) {
	G4cout << " Nuclei fragment: " << G4endl;
	ifrag->printParticle();
      }

      // Copy local G4DynPart to public output 
      aFragment =  new G4DynamicParticle(ifrag->getDynamicParticle());
      theResult.AddSecondary(aFragment); 
    }
  }

  // Report violations of energy, baryon conservation
  if (verboseLevel > 2) {
    if (sumBaryon != 0) {
      G4cout << "ERROR: no baryon number conservation, sum of baryons = "
	     << sumBaryon << G4endl;
    }

    if (sumEnergy > 0.01 ) {
      G4cout << "Kinetic energy conservation violated by "
	     << sumEnergy << " GeV" << G4endl;
    }
     
    G4cout << "Total energy conservation at level ~"
	   << (eInit - eTot) * GeV << " MeV" << G4endl;
    
    if (sumEnergy < -5.0e-5 ) { // 0.05 MeV
      G4cout << "FATAL ERROR: energy created  "
	     << sumEnergy * GeV << " MeV" << G4endl;
    }
  }

  delete bullet;
  delete collider;

  if(target != 0) delete target;
  if(targetH != 0) delete targetH;

  return &theResult;
  }
