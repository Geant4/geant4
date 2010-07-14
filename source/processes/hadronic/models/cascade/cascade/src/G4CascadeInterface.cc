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
// $Id: G4CascadeInterface.cc,v 1.91 2010-07-14 15:41:13 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100414  M. Kelsey -- Check for K0L/K0S before using G4InuclElemPart::type
// 20100418  M. Kelsey -- Reference output particle lists via const-ref, use
//		const_iterator for both.
// 20100428  M. Kelsey -- Use G4InuclParticleNames enum
// 20100429  M. Kelsey -- Change "case gamma:" to "case photon:"
// 20100517  M. Kelsey -- Follow new ctors for G4*Collider family.
// 20100520  M. Kelsey -- Simplify collision loop, move momentum rotations to
//		G4CollisionOutput, copy G4DynamicParticle directly from
//		G4InuclParticle, no switch-block required.
// 20100615  M. Kelsey -- Bug fix: For K0's need ekin in GEANT4 units
// 20100617  M. Kelsey -- Rename "debug_" preprocessor flag to G4CASCADE_DEBUG,
//		and "BERTDEV" to "G4CASCADE_COULOMB_DEV"
// 20100617  M. Kelsey -- Make G4InuclCollider a local data member
// 20100618  M. Kelsey -- Deploy energy-conservation test on final state, with
//		preprocessor flag G4CASCADE_SKIP_ECONS to skip test.
// 20100620  M. Kelsey -- Use new energy-conservation pseudo-collider
// 20100621  M. Kelsey -- Fix compiler warning from GCC 4.5
// 20100624  M. Kelsey -- Fix cascade loop to check nTries every time (had
//		allowed for infinite loop on E-violation); dump event data
//		to output if E-violation exceeds maxTries; use CheckBalance
//		for baryon and charge conservation.
// 20100701  M. Kelsey -- Pass verbosity through to G4CollisionOutput

#include "G4CascadeInterface.hh"
#include "globals.hh"
#include "G4CollisionOutput.hh"
#include "G4DynamicParticle.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4HadronicException.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4KaonZeroShort.hh"
#include "G4KaonZeroLong.hh"
#include "G4LorentzRotation.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4V3DNucleus.hh"
#include <cmath>

using namespace G4InuclParticleNames;

typedef std::vector<G4InuclElementaryParticle>::const_iterator particleIterator;
typedef std::vector<G4InuclNuclei>::const_iterator nucleiIterator;

G4CascadeInterface::G4CascadeInterface(const G4String& nam)
  : G4VIntraNuclearTransportModel(nam), verboseLevel(0) {}


G4CascadeInterface::~G4CascadeInterface() {}
   
G4ReactionProductVector* G4CascadeInterface::Propagate(G4KineticTrackVector* , 
						       G4V3DNucleus* ) {
  return 0;
}

G4HadFinalState* 
G4CascadeInterface::ApplyYourself(const G4HadProjectile& aTrack, 
				  G4Nucleus& theNucleus) {
  if (verboseLevel)
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;

#ifdef G4CASCADE_DEBUG_INTERFACE
  static G4int counter(0);
  counter++;
  G4cerr << "Reaction number "<< counter << " "
	 << aTrack.GetDefinition()->GetParticleName() << " "
	 << aTrack.GetKineticEnergy() << G4endl;
#endif

  theResult.Clear();

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

  G4InuclElementaryParticle* bullet =
    new G4InuclElementaryParticle(momentumBullet, bulletType); 

  // Set target
  G4double theNucleusA = theNucleus.GetN();

  G4InuclParticle* target = 0;
  if (G4int(theNucleusA) == 1)
    target = new G4InuclElementaryParticle(proton);
  else
    target = new G4InuclNuclei(theNucleusA, theNucleus.GetZ());

  if (verboseLevel > 2) {
    G4cout << "Bullet:  " << G4endl;  
    bullet->printParticle();
    G4cout << "Target:  " << G4endl;  
    target->printParticle();
  }

  // Colliders initialisation
  collider.setVerboseLevel(verboseLevel);

  G4CollisionOutput output;
  output.setVerboseLevel(verboseLevel);

  G4CascadeCheckBalance balance(0.05, 0.1, GetModelName());	// Second arg is in GeV
  balance.setVerboseLevel(verboseLevel);

  G4int  maxTries = 100; // maximum tries for inelastic collision to avoid infinite loop
  G4int  nTries   = 0;  // try counter

#ifdef G4CASCADE_COULOMB_DEV
  G4bool coulombOK = false;  // flag for correct Coulomb barrier
#endif

  if (G4int(theNucleusA) == 1) { // special treatment for target H(1,1) (proton)
    G4float cutElastic[32];
    
    cutElastic[proton   ] = 1.0; // GeV
    cutElastic[neutron  ] = 1.0;
    cutElastic[pionPlus ] = 0.6;
    cutElastic[pionMinus] = 0.2;
    cutElastic[pionZero ] = 0.2;
    cutElastic[kaonPlus ] = 0.5;
    cutElastic[kaonMinus] = 0.5;
    cutElastic[kaonZero] = 0.5;
    cutElastic[kaonZeroBar] = 0.5;
    cutElastic[lambda] = 1.0;
    cutElastic[sigmaPlus] = 1.0;
    cutElastic[sigmaZero] = 1.0;
    cutElastic[sigmaMinus] = 1.0;
    cutElastic[xiZero] = 1.0;
    cutElastic[xiMinus] = 1.0;
    
    if (momentumBullet.z() > cutElastic[bulletType]) {
      do {   			// we try to create inelastic interaction
	output.reset();
	collider.collide(bullet, target, output);
	nTries++;
      } while(
	      (nTries < maxTries) &&
	      (output.getOutgoingParticles().size() == 2 &&
	       (output.getOutgoingParticles().begin()->type() == bulletType ||
		output.getOutgoingParticles().begin()->type() == proton)
	       )
	      );
    } else { 		// only elastic collision is energetically possible
      collider.collide(bullet, target, output);
    }
  } else {  			// treat all other targets excepet H(1,1)
    do { 			// we try to create inelastic interaction
      if (verboseLevel > 1)
	G4cout << " Generating cascade attempt " << nTries << G4endl;

      output.reset();
      collider.collide(bullet, target, output);
      nTries++;

      // Check energy conservation; discard result on violation
      balance.collide(bullet, target, output);
      if (verboseLevel > 2) balance.okay();		// Reports violations

#ifdef G4CASCADE_COULOMB_DEV
      coulombOK = false;  		// by default coulomb analysis is OK
      G4double coulumbBarrier = 8.7 * MeV; 
      const std::vector<G4InuclElementaryParticle>& p= output.getOutgoingParticles();
      for (particleIterator ipart = p.begin(); ipart != p.end(); ipart++) {
	if (ipart->type() == proton) {
	  G4double e = ipart->getKineticEnergy()*GeV;
	  
	  // If event with coulomb barrier violation detected -> retry
	  coulombOK |= (e < coulumbBarrier);
	}
      }
#endif
    } while( 
	    (nTries < maxTries) &&  		// conditions for next try
	    ((output.getOutgoingParticles().size()!=0) &&
#ifdef G4CASCADE_COULOMB_DEV
	     (coulombOK) &&
	     ((output.getOutgoingParticles().size() + output.getNucleiFragments().size()) > 2.5)
#else
	     ((output.getOutgoingParticles().size() + output.getNucleiFragments().size()) < 2.5) &&  
	     (output.getOutgoingParticles().begin()->type()==bullet->type())
#endif
	     )
#ifndef G4CASCADE_SKIP_ECONS
	    || (!balance.okay())	// Checks E, p and B conservation
#endif
	     );
  }

  // Check whether repeated attempts have all failed; report and exit
  if (nTries >= maxTries && !balance.okay()) {
    G4cerr << " >>> G4CascadeInterface::ApplyYourself()\n has non-conserving"
	   << " cascade after " << nTries << " attempts." << G4endl;

    G4String throwMsg = "G4CascadeInterface::ApplyYourself() - ";
    if (!balance.energyOkay()) {
      throwMsg += "Energy";
      G4cerr << " Energy conservation violated by " << balance.deltaE()
	     << " GeV (" << balance.relativeE() << ")" << G4endl;
    }

    if (!balance.momentumOkay()) {
      throwMsg += "Momentum";
      G4cerr << " Momentum conservation violated by " << balance.deltaP()
	     << " GeV/c (" << balance.relativeP() << ")" << G4endl;
    }

    if (!balance.baryonOkay()) {
      throwMsg += "Baryon number";
      G4cerr << " Baryon number violated by " << balance.deltaB() << G4endl;
    }

    if (!balance.chargeOkay()) {
      throwMsg += "Charge";
      G4cerr << " Charge conservation violated by " << balance.deltaB()
	     << G4endl;
    }

    G4cout << "\n Final event output, for debugging:"
	   << "\n Bullet:  " << G4endl;  
    bullet->printParticle();
    G4cout << "\n Target:  " << G4endl;  
    target->printParticle();

    output.printCollisionOutput();

    throwMsg += " non-conservation. More info in output.";
    throw G4HadronicException(__FILE__, __LINE__, throwMsg);   // Job ends here!
  }

  // Successful cascade -- clean up and return
  if (verboseLevel > 1) {
    G4cout << " Cascade output after " << nTries << " tries: " << G4endl;
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

      if (!ipart->valid() || ipart->quasi_deutron()) {
        G4cerr << " ERROR: G4CascadeInterface::ApplyYourself incompatible"
	       << " particle type " << outgoingType << G4endl;
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
      if (verboseLevel > 2) {
	G4cout << " Nuclei fragment: " << G4endl;
	ifrag->printParticle();
      }
      
      // Copy local G4DynPart to public output 
      aFragment =  new G4DynamicParticle(ifrag->getDynamicParticle());
      theResult.AddSecondary(aFragment); 
    }
  }

  // Report violations of conservation laws
  balance.collide(bullet, target, output);	// Redo with rotated vectors

  if (verboseLevel > 2) {
    if (!balance.baryonOkay()) {
      G4cerr << "ERROR: no baryon number conservation, sum of baryons = "
             << balance.deltaB() << G4endl;
    }

    if (!balance.chargeOkay()) {
      G4cerr << "ERROR: no charge conservation, sum of charges = "
	     << balance.deltaQ() << G4endl;
    }

    if (std::abs(balance.deltaKE()) > 0.01 ) {	// GeV
      G4cerr << "Kinetic energy conservation violated by "
	     << balance.deltaKE() << " GeV" << G4endl;
    }

    G4double eInit = bullet->getEnergy() + target->getEnergy();
    G4double eFinal = eInit + balance.deltaE();

    G4cout << "Initial energy " << eInit << " final energy " << eFinal
	   << "\nTotal energy conservation at level "
	   << balance.deltaE() * GeV << " MeV" << G4endl;
    
    if (balance.deltaKE() > 5.0e-5 ) { 	// 0.05 MeV
      G4cerr << "FATAL ERROR: kinetic energy created  "
             << balance.deltaKE() * GeV << " MeV" << G4endl;
    }
  }

  delete bullet; bullet=0;
  delete target; target=0;

  return &theResult;
}
