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
// $Id: G4CascadeInterface.cc 71719 2013-06-21 00:01:54Z mkelsey $
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
// 20110224  M. Kelsey -- Add createTarget() for use with Propagate(); split
//		conservation law messages to separate function; begin to add
//		infrastructure code to Propagate.  Move verbose
//		setting to .cc file, and apply to all member objects.
// 20110301  M. Kelsey -- Add copyPreviousCascade() for use with Propagate()  
//		along with new buffers and related particle-conversion  
//		functions.  Encapulate buffer deletion in clear().  Add some
//		diagnostic messages.
// 20110302  M. Kelsey -- Redo diagnostics inside G4CASCADE_DEBUG_INTERFACE
// 20110304  M. Kelsey -- Drop conversion of Propagate() arguments; pass
//		directly to collider for processing.  Rename makeReactionProduct
//		to makeDynamicParticle.
// 20110316  M. Kelsey -- Move kaon-mixing for type-code into G4InuclParticle;
//		add placeholders to capture and use bullet in Propagte
// 20110327  G. Folger -- Set up for E/p checking by G4HadronicProcess in ctor
// 20110328  M. Kelsey -- Modify balance() initialization to match Gunter's
// 20110404  M. Kelsey -- Get primary projectile from base class (ref-03)
// 20110502  M. Kelsey -- Add interface to capture random seeds for user
// 20110719  M. Kelsey -- Use trivialise() in case maximum retries are reached
// 20110720  M. Kelsey -- Discard elastic-cut array (no longer needed),
//		discard local "theFinalState" (in base as "theParticleChange"),
//		Modify createBullet() to set null pointer if bullet unusable,
//		return empty final-state on failures.
//		Fix charge violation report before throwing exception.
// 20110722  M. Kelsey -- In makeDynamicParticle(), allow invalid type codes
//		in order to process, e.g., resonances from Propagate() input
// 20110728  M. Kelsey -- Per V.Ivantchenko, change NoInteraction to return
//		zero particles, but set kinetic energy from projectile.
// 20110801  M. Kelsey -- Make bullet, target point to local buffers, no delete
// 20110802  M. Kelsey -- Use new decay handler for Propagate interface
// 20110922  M. Kelsey -- Follow migration of G4InuclParticle::print(), use
//		G4ExceptionDescription for reporting before throwing exception
// 20120125  M. Kelsey -- In retryInelasticProton() check for empty output
// 20120525  M. Kelsey -- In NoInteraction, check for Ekin<0., set to zero;
//		use SetEnergyChange(0.) explicitly for good final states.
// 20120822  M. Kelsey -- Move envvars to G4CascadeParameters.
// 20130508  D. Wright -- Add support for muon capture
// 20130804  M. Kelsey -- Fix bug #1513 -- "(Z=1)" in boolean expression
// 20140116  M. Kelsey -- Move statics to const data members to avoid weird
//		interactions with MT.
// 20140929  M. Kelsey -- Explicitly call useCascadeDeexcitation() in ctor
// 20150506  M. Kelsey -- Call Initialize() in ctor for master thread only
// 20150608  M. Kelsey -- Label all while loops as terminating.

#include <cmath>
#include <iostream>

#include "G4CascadeInterface.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeCheckBalance.hh"
#include "G4CascadeParameters.hh"
#include "G4CollisionOutput.hh"
#include "G4DecayKineticTracks.hh"
#include "G4DynamicParticle.hh"
#include "G4HadronicException.hh"
#include "G4InuclCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4InuclParticleNames.hh"
#include "G4KaonZeroLong.hh"
#include "G4KaonZeroShort.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Nucleus.hh"
#include "G4ParticleDefinition.hh"
#include "G4ReactionProductVector.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4V3DNucleus.hh"
#include "G4UnboundPN.hh"
#include "G4Dineutron.hh"
#include "G4Diproton.hh"

using namespace G4InuclParticleNames;

typedef std::vector<G4InuclElementaryParticle>::const_iterator particleIterator;
typedef std::vector<G4InuclNuclei>::const_iterator nucleiIterator;


// Constructor and destrutor

G4CascadeInterface::G4CascadeInterface(const G4String& name)
  : G4VIntraNuclearTransportModel(name), 
    randomFile(G4CascadeParameters::randomFile()),
    maximumTries(20), numberOfTries(0),
    collider(new G4InuclCollider), balance(new G4CascadeCheckBalance(name)),
    bullet(0), target(0), output(new G4CollisionOutput) {
  // Set up global objects for master thread or sequential build
  if (G4Threading::IsMasterThread()) Initialize();

  SetEnergyMomentumCheckLevels(5*perCent, 10*MeV);
  balance->setLimits(5*perCent, 10*MeV/GeV);	// Bertini internal units
  this->SetVerboseLevel(G4CascadeParameters::verbose());

  if (G4CascadeParameters::usePreCompound())
    usePreCompoundDeexcitation();
  else
    useCascadeDeexcitation();
}

G4CascadeInterface::~G4CascadeInterface() {
  clear();
  delete collider; collider=0;
  delete balance; balance=0;
  delete output; output=0;
}

void G4CascadeInterface::ModelDescription(std::ostream& outFile) const
{
  outFile << "The Bertini-style cascade implements the inelastic scattering\n"
          << "of hadrons by nuclei.  Nucleons, pions, kaons and hyperons\n"
          << "from 0 to 15 GeV may be used as projectiles in this model.\n"
          << "Final state hadrons are produced by a classical cascade\n"
          << "consisting of individual hadron-nucleon scatterings which use\n"
          << "free-space partial cross sections, corrected for various\n"
          << "nuclear medium effects.  The target nucleus is modeled as a\n"
          << "set of 1, 3 or 6 spherical shells, in which scattered hadrons\n"
          << "travel in straight lines until they are reflected from or\n"
          << "transmitted through shell boundaries.\n";
}

void G4CascadeInterface::DumpConfiguration(std::ostream& outFile) const {
  G4CascadeParameters::DumpConfiguration(outFile);
}

void G4CascadeInterface::clear() {
  bullet=0;
  target=0;
}


// Initialize shared objects for use across multiple threads

void G4CascadeInterface::Initialize() {
  G4ParticleDefinition* pn = G4UnboundPN::Definition();
  G4ParticleDefinition* nn = G4Dineutron::Definition();
  G4ParticleDefinition* pp = G4Diproton::Definition();
  const G4CascadeChannel* ch = G4CascadeChannelTables::GetTable(0);
  if (!ch || !pn || !nn || !pp) return;		// Avoid "unused variables"
}


// Select post-cascade processing (default will be CascadeDeexcitation)
// NOTE:  Currently just calls through to Collider, in future will do something

void G4CascadeInterface::useCascadeDeexcitation() {
  collider->useCascadeDeexcitation();
}

void G4CascadeInterface::usePreCompoundDeexcitation() {
  collider->usePreCompoundDeexcitation();
}


// Apply verbosity to all member objects (overrides base class version)

void G4CascadeInterface::SetVerboseLevel(G4int verbose) {
  G4HadronicInteraction::SetVerboseLevel(verbose);
  collider->setVerboseLevel(verbose);
  balance->setVerboseLevel(verbose);
  output->setVerboseLevel(verbose);
}


// Test whether inputs are valid for this model

G4bool G4CascadeInterface::IsApplicable(const G4HadProjectile& aTrack,
					G4Nucleus& /* theNucleus */) {
  return IsApplicable(aTrack.GetDefinition());  
}

G4bool G4CascadeInterface::IsApplicable(const G4ParticleDefinition* aPD) const {
  if (aPD->GetAtomicMass() > 1) return true;		// Nuclei are okay

  // Valid particle and have interactions available
  G4int type = G4InuclElementaryParticle::type(aPD);
  return (G4CascadeChannelTables::GetTable(type));
}


// Main Actions

G4HadFinalState* 
G4CascadeInterface::ApplyYourself(const G4HadProjectile& aTrack, 
				  G4Nucleus& theNucleus) {
  if (verboseLevel)
    G4cout << " >>> G4CascadeInterface::ApplyYourself" << G4endl;

  if (aTrack.GetKineticEnergy() < 0.) {
    G4cerr << " >>> G4CascadeInterface got negative-energy track: "
	   << aTrack.GetDefinition()->GetParticleName() << " Ekin = "
	   << aTrack.GetKineticEnergy() << G4endl;
  }

#ifdef G4CASCADE_DEBUG_INTERFACE
  static G4int counter(0);
  counter++;
  G4cerr << "Reaction number "<< counter << " "
	 << aTrack.GetDefinition()->GetParticleName() << " "
	 << aTrack.GetKineticEnergy() << " MeV" << G4endl;
#endif

  if (!randomFile.empty()) {		// User requested random-seed capture
    if (verboseLevel>1) 
      G4cout << " Saving random engine state to " << randomFile << G4endl;
    CLHEP::HepRandom::saveEngineStatus(randomFile);
  }

  theParticleChange.Clear();
  clear();

  // Abort processing if no interaction is possible
  if (!IsApplicable(aTrack, theNucleus)) {
    if (verboseLevel) G4cerr << " No interaction possible " << G4endl;
    return NoInteraction(aTrack, theNucleus);
  }

  // Make conversion between native Geant4 and Bertini cascade classes.
  if (!createBullet(aTrack)) {
    if (verboseLevel) G4cerr << " Unable to create usable bullet" << G4endl;
    return NoInteraction(aTrack, theNucleus);
  }

  if (!createTarget(theNucleus)) {
    if (verboseLevel) G4cerr << " Unable to create usable target" << G4endl;
    return NoInteraction(aTrack, theNucleus);
  }

  // Different retry conditions for proton target vs. nucleus
  const G4bool isHydrogen = (theNucleus.GetA_asInt() == 1);

  numberOfTries = 0;
  do {   			// we try to create inelastic interaction
    if (verboseLevel > 1)
      G4cout << " Generating cascade attempt " << numberOfTries << G4endl;
    
    output->reset();
    collider->collide(bullet, target, *output);
    balance->collide(bullet, target, *output);
    
    numberOfTries++;
    /* Loop checking 08.06.2015 MHK */
  } while ( isHydrogen ? retryInelasticProton() : retryInelasticNucleus() );

  // Null event if unsuccessful
  if (numberOfTries >= maximumTries) {
    if (verboseLevel) 
      G4cout << " Cascade aborted after trials " << numberOfTries << G4endl;
    return NoInteraction(aTrack, theNucleus);
  }

  // Abort job if energy or momentum are not conserved
  if (!balance->okay()) {
    throwNonConservationFailure();
    return NoInteraction(aTrack, theNucleus);
  }

  // Successful cascade -- clean up and return
  if (verboseLevel) {
    G4cout << " Cascade output after trials " << numberOfTries << G4endl;
    if (verboseLevel > 1) output->printCollisionOutput();
  }

  // Rotate event to put Z axis along original projectile direction
  // Removed by DHW to fix bug #1990
  // output->rotateEvent(bulletInLabFrame);

  copyOutputToHadronicResult();

  // Report violations of conservation laws in original frame
  checkFinalResult();

  // Clean up and return final result;
  clear();
/*
  G4int nSec = theParticleChange.GetNumberOfSecondaries();
  for (G4int i = 0; i < nSec; i++) {
    G4HadSecondary* sec = theParticleChange.GetSecondary(i);
    G4DynamicParticle* dp = sec->GetParticle();
    if (dp->GetDefinition()->GetParticleName() == "neutron") 
      G4cout << dp->GetDefinition()->GetParticleName() << " has "
             << dp->GetKineticEnergy()/MeV << " MeV " << G4endl;
  }
*/
  return &theParticleChange;
}

G4ReactionProductVector* 
G4CascadeInterface::Propagate(G4KineticTrackVector* theSecondaries,
			      G4V3DNucleus* theNucleus) {
  if (verboseLevel) G4cout << " >>> G4CascadeInterface::Propagate" << G4endl;

#ifdef G4CASCADE_DEBUG_INTERFACE
  if (verboseLevel>1) {
    G4cout << " G4V3DNucleus A " << theNucleus->GetMassNumber()
	   << " Z " << theNucleus->GetCharge()
	   << "\n " << theSecondaries->size() << " secondaries:" << G4endl;
    for (size_t i=0; i<theSecondaries->size(); i++) {
      G4KineticTrack* kt = (*theSecondaries)[i];
      G4cout << " " << i << ": " << kt->GetDefinition()->GetParticleName() 
	     << " p " << kt->Get4Momentum() << " @ " << kt->GetPosition()
	     << " t " << kt->GetFormationTime() << G4endl;
    }
  }
#endif

  if (!randomFile.empty()) {		// User requested random-seed capture
    if (verboseLevel>1) 
      G4cout << " Saving random engine state to " << randomFile << G4endl;
    CLHEP::HepRandom::saveEngineStatus(randomFile);
  }

  theParticleChange.Clear();
  clear();

  // Process input secondaries list to eliminate resonances
  G4DecayKineticTracks decay(theSecondaries);

  // NOTE:  Requires 9.4-ref-03 mods to base class and G4TheoFSGenerator
  const G4HadProjectile* projectile = GetPrimaryProjectile();
  if (projectile) createBullet(*projectile);

  if (!createTarget(theNucleus)) {
    if (verboseLevel) G4cerr << " Unable to create usable target" << G4endl;
    return 0;	// FIXME:  This will cause a segfault later
  }

  numberOfTries = 0;
  do {
    if (verboseLevel > 1)
      G4cout << " Generating rescatter attempt " << numberOfTries << G4endl;

    output->reset();
    collider->rescatter(bullet, theSecondaries, theNucleus, *output);
    balance->collide(bullet, target, *output);

    numberOfTries++;
    // FIXME:  retry checks will SEGFAULT until we can define the bullet!
  } while (retryInelasticNucleus());	/* Loop checking 08.06.2015 MHK */

  // Check whether repeated attempts have all failed; report and exit
  if (numberOfTries >= maximumTries && !balance->okay()) {
    throwNonConservationFailure();	// This terminates the job
  }

  // Successful cascade -- clean up and return
  if (verboseLevel) {
    G4cout << " Cascade rescatter after trials " << numberOfTries << G4endl;
    if (verboseLevel > 1) output->printCollisionOutput();
  }

  // Does calling code take ownership?  I hope so!
  G4ReactionProductVector* propResult = copyOutputToReactionProducts();

  // Clean up and and return final result
  clear();
  return propResult;
}


// Replicate input particles onto output

G4HadFinalState* 
G4CascadeInterface::NoInteraction(const G4HadProjectile& aTrack,
				  G4Nucleus& /*theNucleus*/) {
  if (verboseLevel) 
    G4cout << " >>> G4CascadeInterface::NoInteraction" << G4endl;

  theParticleChange.Clear();
  theParticleChange.SetStatusChange(isAlive);

  G4double ekin = aTrack.GetKineticEnergy()>0. ? aTrack.GetKineticEnergy() : 0.;
  theParticleChange.SetEnergyChange(ekin);	// Protect against rounding

  return &theParticleChange;
}


// Convert input projectile to Bertini internal object

G4bool G4CascadeInterface::createBullet(const G4HadProjectile& aTrack) {
  const G4ParticleDefinition* trkDef = aTrack.GetDefinition();
  G4int bulletType = 0;			// For elementary particles
  G4int bulletA = 0, bulletZ = 0;	// For nucleus projectile

  if (trkDef->GetAtomicMass() <= 1) {
    bulletType = G4InuclElementaryParticle::type(trkDef);
  } else {
    bulletA = trkDef->GetAtomicMass();
    bulletZ = trkDef->GetAtomicNumber();
  }

  if (0 == bulletType && 0 == bulletA*bulletZ) {
    if (verboseLevel) {
      G4cerr << " G4CascadeInterface: " << trkDef->GetParticleName()
	     << " not usable as bullet." << G4endl;
    }
    bullet = 0;
    return false;
  }

  // Code momentum and energy -- Bertini wants z-axis and GeV units
  G4LorentzVector projectileMomentum = aTrack.Get4Momentum()/GeV;

  // Rotation/boost to get from z-axis back to original frame
  // According to bug report #1990 this rotation is unnecessary and causes
  // irreproducibility.  Verifed and fixed by DHW 27 Nov 2017
  // bulletInLabFrame = G4LorentzRotation::IDENTITY;	// Initialize
  // bulletInLabFrame.rotateZ(-projectileMomentum.phi());
  // bulletInLabFrame.rotateY(-projectileMomentum.theta());
  // bulletInLabFrame.invert();

  G4LorentzVector momentumBullet(0., 0., projectileMomentum.rho(),
				 projectileMomentum.e());
  
  if (G4InuclElementaryParticle::valid(bulletType)) {
    hadronBullet.fill(momentumBullet, bulletType);
    bullet = &hadronBullet;
  } else {
    nucleusBullet.fill(momentumBullet, bulletA, bulletZ);
    bullet = &nucleusBullet;
  }

  if (verboseLevel > 2) G4cout << "Bullet:  \n" << *bullet << G4endl;  

  return true;
}


// Convert input nuclear target to Bertini internal object

G4bool G4CascadeInterface::createTarget(G4Nucleus& theNucleus) {
  return createTarget(theNucleus.GetA_asInt(), theNucleus.GetZ_asInt());
}

G4bool G4CascadeInterface::createTarget(G4V3DNucleus* theNucleus) {
  return createTarget(theNucleus->GetMassNumber(), theNucleus->GetCharge());
}

G4bool G4CascadeInterface::createTarget(G4int A, G4int Z) {
  if (A > 1) {
    nucleusTarget.fill(A, Z);
    target = &nucleusTarget;
  } else {
    hadronTarget.fill(0., (Z==1?proton:neutron));
    target = &hadronTarget;
  }

  if (verboseLevel > 2) G4cout << "Target:  \n" << *target << G4endl;

  return true;		// Right now, target never fails
}


// Convert Bertini particle to output (G4DynamicParticle)

G4DynamicParticle* G4CascadeInterface::
makeDynamicParticle(const G4InuclElementaryParticle& iep) const {
  G4int outgoingType = iep.type();
  
  if (iep.quasi_deutron()) {
    G4cerr << " ERROR: G4CascadeInterface incompatible particle type "
	   << outgoingType << G4endl;
    return 0;
  }
  
  // Copy local G4DynPart to public output (handle kaon mixing specially)
  if (outgoingType == kaonZero || outgoingType == kaonZeroBar) {
    G4ThreeVector momDir = iep.getMomentum().vect().unit();
    G4double ekin = iep.getKineticEnergy()*GeV;	// Bertini -> G4 units
    
    G4ParticleDefinition* pd = G4KaonZeroShort::Definition();
    if (G4UniformRand() > 0.5) pd = G4KaonZeroLong::Definition();
    
    return new G4DynamicParticle(pd, momDir, ekin);
  } else {
    return new G4DynamicParticle(iep.getDynamicParticle());
  }
  
  return 0;	// Should never get here!
}
 
G4DynamicParticle* 
G4CascadeInterface::makeDynamicParticle(const G4InuclNuclei& inuc) const {
  if (verboseLevel > 2) {
    G4cout << " Nuclei fragment: \n" << inuc << G4endl;
  }
  
  // Copy local G4DynPart to public output 
  return new G4DynamicParticle(inuc.getDynamicParticle());
}


// Transfer Bertini internal final state to hadronics interface

void G4CascadeInterface::copyOutputToHadronicResult() {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeInterface::copyOutputToHadronicResult" << G4endl;

  const std::vector<G4InuclNuclei>& outgoingNuclei = output->getOutgoingNuclei();
  const std::vector<G4InuclElementaryParticle>& particles = output->getOutgoingParticles();

  theParticleChange.SetStatusChange(stopAndKill);
  theParticleChange.SetEnergyChange(0.);

  // Get outcoming particles
  if (!particles.empty()) { 
    particleIterator ipart = particles.begin();
    for (; ipart != particles.end(); ipart++) {
      theParticleChange.AddSecondary(makeDynamicParticle(*ipart));
    }
  }

  // get nuclei fragments
  if (!outgoingNuclei.empty()) { 
    nucleiIterator ifrag = outgoingNuclei.begin();
    for (; ifrag != outgoingNuclei.end(); ifrag++) {
      theParticleChange.AddSecondary(makeDynamicParticle(*ifrag)); 
    }
  }
}

G4ReactionProductVector* G4CascadeInterface::copyOutputToReactionProducts() {
  if (verboseLevel > 1)
    G4cout << " >>> G4CascadeInterface::copyOutputToReactionProducts" << G4endl;

  const std::vector<G4InuclElementaryParticle>& particles = output->getOutgoingParticles();
  const std::vector<G4InuclNuclei>& fragments = output->getOutgoingNuclei();

  G4ReactionProductVector* propResult = new G4ReactionProductVector;

  G4ReactionProduct* rp = 0;	// Buffers to create outgoing tracks
  G4DynamicParticle* dp = 0;

  // Get outcoming particles
  if (!particles.empty()) { 
    particleIterator ipart = particles.begin();
    for (; ipart != particles.end(); ipart++) {
      rp = new G4ReactionProduct;
      dp = makeDynamicParticle(*ipart);
      (*rp) = (*dp);		// This does all the necessary copying
      propResult->push_back(rp);
      delete dp;
    }
  }

  // get nuclei fragments
  if (!fragments.empty()) { 
    nucleiIterator ifrag = fragments.begin();
    for (; ifrag != fragments.end(); ifrag++) {
      rp = new G4ReactionProduct;
      dp = makeDynamicParticle(*ifrag);
      (*rp) = (*dp);		// This does all the necessary copying
      propResult->push_back(rp);
      delete dp;
    }
  }

  return propResult;
}


// Report violations of conservation laws in original frame

void G4CascadeInterface::checkFinalResult() {
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

#ifdef G4CASCADE_DEBUG_INTERFACE
  // Report on all retry conditions, in order of return logic
  G4cout << " retryInelasticProton: number of Tries "
	 << ((numberOfTries < maximumTries) ? "RETRY (t)" : "EXIT (f)")
	 << "\n retryInelasticProton: AND collision type ";
  if (out.empty()) G4cout << "FAILED" << G4endl;
  else {
    G4cout << (out.size() == 2 ? "ELASTIC (t)" : "INELASTIC (f)")
	   << "\n retryInelasticProton: AND Leading particles bullet "
	   << (out.size() >= 2 &&
	       (out[0].getDefinition() == bullet->getDefinition() ||
		out[1].getDefinition() == bullet->getDefinition())
	       ? "YES (t)" : "NO (f)")
	   << G4endl;
  }
#endif

  return ( (numberOfTries < maximumTries) &&
	   (out.empty() ||
	    (out.size() == 2 &&
	     (out[0].getDefinition() == bullet->getDefinition() ||
	      out[1].getDefinition() == bullet->getDefinition())))
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

#ifdef G4CASCADE_DEBUG_INTERFACE
  // Report on all retry conditions, in order of return logic
  G4cout << " retryInelasticNucleus: numberOfTries "
	 << ((numberOfTries < maximumTries) ? "RETRY (t)" : "EXIT (f)")
	 << "\n retryInelasticNucleus: AND outputParticles "
	 << ((npart != 0) ? "NON-ZERO (t)" : "EMPTY (f)")
#ifdef G4CASCADE_COULOMB_DEV
	 << "\n retryInelasticNucleus: AND coulombBarrier (COULOMB_DEV) "
	 << (coulombBarrierViolation() ? "VIOLATED (t)" : "PASSED (f)")
	 << "\n retryInelasticNucleus: AND collision type (COULOMB_DEV) "
	 << ((npart+nfrag > 2) ? "INELASTIC (t)" : "ELASTIC (f)")
#else
	 << "\n retryInelasticNucleus: AND collsion type "
	 << ((npart+nfrag < 3) ? "ELASTIC (t)" : "INELASTIC (f)")
	 << "\n retryInelasticNucleus: AND Leading particle bullet "
	 << ((firstOut == bullet->getDefinition()) ? "YES (t)" : "NO (f)")
#endif
	 << "\n retryInelasticNucleus: OR conservation "
	 << (!balance->okay() ? "FAILED (t)" : "PASSED (f)")
	 << G4endl;
#endif

  return ( (numberOfTries < maximumTries) &&
	   ( ((npart != 0) &&
#ifdef G4CASCADE_COULOMB_DEV
	      (coulombBarrierViolation() && npart+nfrag > 2)
#else
	      (npart+nfrag < 3 && firstOut == bullet->getDefinition())
#endif
             )
#ifndef G4CASCADE_SKIP_ECONS
	     || (!balance->okay())
#endif
           )
	 );
}


// Terminate job in case of persistent non-conservation
// FIXME:  Need to migrate to G4ExceptionDescription

void G4CascadeInterface::throwNonConservationFailure() {
  // NOTE:  Once G4HadronicException is changed, use the following line!
  // G4ExceptionDescription errInfo;
  std::ostream& errInfo = G4cerr;

  errInfo << " >>> G4CascadeInterface has non-conserving"
	  << " cascade after " << numberOfTries << " attempts." << G4endl;

  G4String throwMsg = "G4CascadeInterface - ";
  if (!balance->energyOkay()) {
    throwMsg += "Energy";
    errInfo << " Energy conservation violated by " << balance->deltaE()
	    << " GeV (" << balance->relativeE() << ")" << G4endl;
  }
  
  if (!balance->momentumOkay()) {
    throwMsg += "Momentum";
    errInfo << " Momentum conservation violated by " << balance->deltaP()
	    << " GeV/c (" << balance->relativeP() << ")" << G4endl;
  }
  
  if (!balance->baryonOkay()) {
    throwMsg += "Baryon number";
    errInfo << " Baryon number violated by " << balance->deltaB() << G4endl;
  }
  
  if (!balance->chargeOkay()) {
    throwMsg += "Charge";
    errInfo << " Charge conservation violated by " << balance->deltaQ()
	    << G4endl;
  }

  errInfo << " Final event output, for debugging:\n"
	 << " Bullet:  \n" << *bullet << G4endl
	 << " Target:  \n" << *target << G4endl;
  output->printCollisionOutput(errInfo);
  
  throwMsg += " non-conservation. More info in output.";
  throw G4HadronicException(__FILE__, __LINE__, throwMsg);   // Job ends here!
}
