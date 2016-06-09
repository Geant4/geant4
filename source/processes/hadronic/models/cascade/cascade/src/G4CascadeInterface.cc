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
//
// $Id: G4CascadeInterface.cc,v 1.105 2010-12-15 07:40:15 gunter Exp $
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
// 20100714  M. Kelsey -- Report number of iterations before success
// 20100720  M. Kelsey -- Use G4CASCADE_SKIP_ECONS flag for reporting
// 20100723  M. Kelsey -- Move G4CollisionOutput to .hh file for reuse
// 20100914  M. Kelsey -- Migrate to integer A and Z
// 20100916  M. Kelsey -- Simplify ApplyYourself() by encapsulating code blocks
//		into numerous functions; make data-member colliders pointers;
//		provide support for projectile nucleus
// 20100919  M. Kelsey -- Fix incorrect logic in retryInelasticNucleus()
// 20100922  M. Kelsey -- Add functions to select de-excitation method
// 20100924  M. Kelsey -- Migrate to "OutgoingNuclei" names in CollisionOutput 
// 20111117  J. Apostolakis -- Reduce use of Fatal Exception for bad E/p on
//		Hydrogen 

#include "G4CascadeInterface.hh"
#include "globals.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CollisionOutput.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronicException.hh"
#include "G4InuclCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4Track.hh"
#include "G4V3DNucleus.hh"
#include <cmath>

using namespace G4InuclParticleNames;

typedef std::vector<G4InuclElementaryParticle>::const_iterator particleIterator;
typedef std::vector<G4InuclNuclei>::const_iterator nucleiIterator;


// Maximum number of iterations allowed for inelastic collision attempts

const G4int G4CascadeInterface::maximumTries = 100;


// Constructor and destrutor

G4CascadeInterface::G4CascadeInterface(const G4String& name)
  : G4VIntraNuclearTransportModel(name),
    verboseLevel(0), numberOfTries(0),
    collider(new G4InuclCollider), 
    balance(new G4CascadeCheckBalance(0.05, 0.1, name)),
    bullet(0), target(0), output(new G4CollisionOutput) {
  initializeElasticCuts();
}

G4CascadeInterface::~G4CascadeInterface() {
  delete collider; collider=0;
  delete balance; balance=0;
  delete bullet; bullet=0;
  delete target; target=0;
  delete output; output=0;
}

// Fill sparse array with minimum momenta for inelastic on hydrogen

void G4CascadeInterface::initializeElasticCuts() {
  cutElastic[proton   ]   = 1.0;	// Bertini uses GeV for everything
  cutElastic[neutron  ]   = 1.0;
  cutElastic[pionPlus ]   = 0.6;
  cutElastic[pionMinus]   = 0.2;
  cutElastic[pionZero ]   = 0.2;
  cutElastic[kaonPlus ]   = 0.5;
  cutElastic[kaonMinus]   = 0.5;
  cutElastic[kaonZero]    = 0.5;
  cutElastic[kaonZeroBar] = 0.5;
  cutElastic[lambda]      = 1.0;
  cutElastic[sigmaPlus]   = 1.0;
  cutElastic[sigmaZero]   = 1.0;
  cutElastic[sigmaMinus]  = 1.0;
  cutElastic[xiZero]      = 1.0;
  cutElastic[xiMinus]     = 1.0;
}


// Select post-cascade processing (default will be CascadeDeexcitation)
// NOTE:  Currently just calls through to Collider, in future will do something

void G4CascadeInterface::useCascadeDeexcitation() {
  collider->useCascadeDeexcitation();
}

void G4CascadeInterface::usePreCompoundDeexcitation() {
  collider->usePreCompoundDeexcitation();
}

// Main Actions

G4ReactionProductVector* 
G4CascadeInterface::Propagate(G4KineticTrackVector*, G4V3DNucleus* ) {
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

  createBullet(aTrack);
  createTarget(theNucleus);

  if (verboseLevel > 2) {
    G4cout << "Bullet:  " << G4endl;  
    bullet->printParticle();
    G4cout << "Target:  " << G4endl;  
    target->printParticle();
  }

  // Colliders initialisation
  collider->setVerboseLevel(verboseLevel);
  balance->setVerboseLevel(verboseLevel);
  output->setVerboseLevel(verboseLevel);

  numberOfTries = 0;

  // special treatment for target H(1,1) (proton)
  if (theNucleus.GetA_asInt() == 1) {
    G4int btype = dynamic_cast<G4InuclElementaryParticle*>(bullet)->type();

    if (bullet->getMomModule() <= cutElastic[btype]) {
      collider->collide(bullet, target, *output);	// Only elastic
    } else {
      do {   			// we try to create inelastic interaction
	if (verboseLevel > 1)
	  G4cout << " Generating cascade attempt " << numberOfTries << G4endl;

	output->reset();
	collider->collide(bullet, target, *output);
	balance->collide(bullet, target, *output);

	numberOfTries++;
      } while(retryInelasticProton());
    }
    // Check whether repeated attempts have all failed; 
    if (numberOfTries >= maximumTries && !balance->okay()) {
      // If E/p problem only report - this is a known issue for Hydrogen

      if(!(balance->momentumOkay() || balance->energyOkay()) 
	 && balance->baryonOkay() && balance->chargeOkay()) {
	const  G4double minDiffEP= 0.001;  // GeV or GeV/c
	static G4double EmaxDiffFound = -100.00; // 
	static G4double PmaxDiffFound = -100.00; // 
        static int      countBadEP=0;
	G4bool badEnergy=   (!balance->energyOkay()) && (balance->deltaE() > minDiffEP);
	G4bool badMomentum= (!balance->momentumOkay()) && (balance->deltaP() > minDiffEP);

        countBadEP++;
	if (badEnergy || badMomentum) {
	  G4cerr << "Warning: " << __FILE__ << ":" 
		 << " In Bertini model - problem with conservation law for H target. " 
		 << G4endl;
	  if (badEnergy) {
	    G4double diffE= balance->deltaE();
	    if ((diffE > EmaxDiffFound) || ((countBadEP%100)==1)) {
	      G4cerr << " Energy conservation violated by " << balance->deltaE()
		     << " GeV (" << balance->relativeE() << ")" << G4endl;
	    }
	    EmaxDiffFound = std::max( diffE, EmaxDiffFound); 
	  }
	  if (badMomentum) {
	    G4double diffP= balance->deltaP();
	    if ((diffP > PmaxDiffFound) || ((countBadEP%100)==1)) {	  
	      G4cerr << " Momentum conservation violated by " << balance->deltaP()
		     << " GeV/c (" << balance->relativeP() << ")" << G4endl;
	    }
	    PmaxDiffFound = std::max( diffP, EmaxDiffFound); 
	  }
	}
      } else
	throwNonConservationFailure();	// This terminates the job
    }
  } else {  			// treat all other targets excepet H(1,1)
    do { 			// we try to create inelastic interaction
      if (verboseLevel > 1)
	G4cout << " Generating cascade attempt " << numberOfTries << G4endl;

      output->reset();
      collider->collide(bullet, target, *output);
      balance->collide(bullet, target, *output);

      numberOfTries++;
    } while (retryInelasticNucleus());

    // Check whether repeated attempts have all failed; report and exit
    if (numberOfTries >= maximumTries && !balance->okay()) {
      throwNonConservationFailure();	// This terminates the job
    }
  }

  // Successful cascade -- clean up and return
  if (verboseLevel) {
    G4cout << " Cascade output after trials " << numberOfTries << G4endl;
    if (verboseLevel > 1) output->printCollisionOutput();
  }

  // Rotate event to put Z axis along original projectile direction
  output->rotateEvent(bulletInLabFrame);

  copyOutputToHadronicResult();

  // Report violations of conservation laws in original frame
  balance->collide(bullet, target, *output);

  if (verboseLevel > 2) {
    if (!balance->baryonOkay()) {
      G4cerr << "ERROR: no baryon number conservation, sum of baryons = "
             << balance->deltaB() << G4endl;
    }

    if (!balance->chargeOkay()) {
      G4cerr << "ERROR: no charge conservation, sum of charges = "
	     << balance->deltaQ() << G4endl;
    }

    if (std::abs(balance->deltaKE()) > 0.01 ) {	// GeV
      G4cerr << "Kinetic energy conservation violated by "
	     << balance->deltaKE() << " GeV" << G4endl;
    }

    G4double eInit = bullet->getEnergy() + target->getEnergy();
    G4double eFinal = eInit + balance->deltaE();

    G4cout << "Initial energy " << eInit << " final energy " << eFinal
	   << "\nTotal energy conservation at level "
	   << balance->deltaE() * GeV << " MeV" << G4endl;
    
    if (balance->deltaKE() > 5.0e-5 ) { 	// 0.05 MeV
      G4cerr << "FATAL ERROR: kinetic energy created  "
             << balance->deltaKE() * GeV << " MeV" << G4endl;
    }
  }

  delete bullet; bullet=0;
  delete target; target=0;

  return &theResult;
}


// Convert input projectile to Bertini internal object

void G4CascadeInterface::createBullet(const G4HadProjectile& aTrack) {
  const G4ParticleDefinition* trkDef = aTrack.GetDefinition();

  G4int bulletType = 0;			// For elementary particles
  G4int bulletA = 0, bulletZ = 0;	// For nucleus projectile

  if (trkDef->GetAtomicMass() <= 1) {
    if (trkDef == G4KaonZeroLong::KaonZeroLong() ||
	trkDef == G4KaonZeroShort::KaonZeroShort() )
      bulletType = (G4UniformRand() > 0.5) ? kaonZero : kaonZeroBar;
    else 
      bulletType = G4InuclElementaryParticle::type(trkDef);
  } else {
    bulletA = trkDef->GetAtomicMass();
    bulletZ = trkDef->GetAtomicNumber();
  }

  // Code momentum and energy -- Bertini wants z-axis and GeV units
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum()/GeV;
  
  // Rrotation/boost to get from z-axis back to original frame
  bulletInLabFrame = G4LorentzRotation::IDENTITY;	// Initialize
  bulletInLabFrame.rotateZ(-projectileMomentum.phi());
  bulletInLabFrame.rotateY(-projectileMomentum.theta());
  bulletInLabFrame.invert();
  
  G4LorentzVector momentumBullet(0., 0., projectileMomentum.rho(),
				 projectileMomentum.e());

  if (bulletType > 0)
    bullet = new G4InuclElementaryParticle(momentumBullet, bulletType);
  else
    bullet = new G4InuclNuclei(momentumBullet, bulletA, bulletZ);

  return;
}


// Convert input nuclear target to Bertini internal object

void G4CascadeInterface::createTarget(G4Nucleus& theNucleus) {
  G4int theNucleusA = theNucleus.GetA_asInt();
  G4int theNucleusZ = theNucleus.GetZ_asInt();

  if (theNucleusA == 1)
    target = new G4InuclElementaryParticle((theNucleusZ==1)?proton:neutron);
  else
    target = new G4InuclNuclei(theNucleusA, theNucleusZ);

  return;
}


// Transfer Bertini internal final state to hadronics interface

void G4CascadeInterface::copyOutputToHadronicResult() {
  const std::vector<G4InuclNuclei>& outgoingNuclei = output->getOutgoingNuclei();
  const std::vector<G4InuclElementaryParticle>& particles = output->getOutgoingParticles();

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
  if (!outgoingNuclei.empty()) { 
    nucleiIterator ifrag = outgoingNuclei.begin();
    for (; ifrag != outgoingNuclei.end(); ifrag++) {
      if (verboseLevel > 2) {
	G4cout << " Nuclei fragment: " << G4endl;
	ifrag->printParticle();
      }
      
      // Copy local G4DynPart to public output 
      aFragment =  new G4DynamicParticle(ifrag->getDynamicParticle());
      theResult.AddSecondary(aFragment); 
    }
  }
}


// Evaluate whether any outgoing particles penetrated Coulomb barrier

G4bool G4CascadeInterface::coulombBarrierViolation() const {
  G4bool violated = false;  		// by default coulomb analysis is OK

  const G4double coulumbBarrier = 8.7 * MeV/GeV; 	// Bertini uses GeV

  const std::vector<G4InuclElementaryParticle>& p =
    output->getOutgoingParticles();

  for (particleIterator ipart=p.begin(); ipart != p.end(); ipart++) {
    if (ipart->type() == proton) {
      violated |= (ipart->getKineticEnergy() < coulumbBarrier);
    }
  }

  return violated;
}

// Check whether inelastic collision on proton failed

G4bool G4CascadeInterface::retryInelasticProton() const {
  const std::vector<G4InuclElementaryParticle>& out =
    output->getOutgoingParticles();

  return ( (numberOfTries < maximumTries) &&
	   (out.size() == 2) &&
	   (out[0].getDefinition() == bullet->getDefinition() ||
	    out[1].getDefinition() == bullet->getDefinition())
	   );
}

// Check whether generic inelastic collision failed
// NOTE:  some conditions are set by compiler flags

G4bool G4CascadeInterface::retryInelasticNucleus() const {
  // Quantities necessary for conditional block below
  G4int npart = output->numberOfOutgoingParticles();
  G4int nfrag = output->numberOfOutgoingNuclei();

  const G4ParticleDefinition* firstOut = (npart == 0) ? 0 :
    output->getOutgoingParticles().begin()->getDefinition();

  return ( ((numberOfTries < maximumTries) &&
	    (npart != 0) &&
#ifdef G4CASCADE_COULOMB_DEV
	    (coulombBarrierViolation() && npart+nfrag > 2)
#else
	    (npart+nfrag < 3 && firstOut == bullet->getDefinition())
#endif
	    )
#ifndef G4CASCADE_SKIP_ECONS
	   || (!balance->okay())
#endif
	   );
}


// Terminate job in case of persistent non-conservation

void G4CascadeInterface::throwNonConservationFailure() {
  G4cerr << " >>> G4CascadeInterface::ApplyYourself()\n has non-conserving"
	 << " cascade after " << numberOfTries << " attempts." << G4endl;

  G4String throwMsg = "G4CascadeInterface::ApplyYourself() - ";
  if (!balance->energyOkay()) {
    throwMsg += "Energy";
    G4cerr << " Energy conservation violated by " << balance->deltaE()
	   << " GeV (" << balance->relativeE() << ")" << G4endl;
  }
  
  if (!balance->momentumOkay()) {
    throwMsg += "Momentum";
    G4cerr << " Momentum conservation violated by " << balance->deltaP()
	   << " GeV/c (" << balance->relativeP() << ")" << G4endl;
  }
  
  if (!balance->baryonOkay()) {
    throwMsg += "Baryon number";
    G4cerr << " Baryon number violated by " << balance->deltaB() << G4endl;
  }
  
  if (!balance->chargeOkay()) {
    throwMsg += "Charge";
    G4cerr << " Charge conservation violated by " << balance->deltaB()
	   << G4endl;
  }
  
  G4cout << "\n Final event output, for debugging:"
	 << "\n Bullet:  " << G4endl;  
  bullet->printParticle();
  G4cout << "\n Target:  " << G4endl;  
  target->printParticle();
  
  output->printCollisionOutput();
  
  throwMsg += " non-conservation. More info in output.";
  throw G4HadronicException(__FILE__, __LINE__, throwMsg);   // Job ends here!
}
