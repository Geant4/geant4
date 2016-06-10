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
// $Id: G4IntraNucleiCascader.cc 67746 2013-03-05 21:11:14Z mkelsey $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100307  M. Kelsey -- Bug fix: momentum_out[0] should be momentum_out.e()
// 20100309  M. Kelsey -- Eliminate some unnecessary std::pow()
// 20100407  M. Kelsey -- Pass "all_particles" as argument to initializeCascad,
//		following recent change to G4NucleiModel.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members
// 20100616  M. Kelsey -- Add reporting of final residual particle
// 20100617  M. Kelsey -- Remove "RUN" preprocessor flag and all "#else" code,
//		pass verbosity to collider.  Make G4NucleiModel a data member,
//		instead of creating and deleting on every cycle.
// 20100620  M. Kelsey -- Improved diagnostic messages.  Simplify kinematics
//		of recoil nucleus.
// 20100622  M. Kelsey -- Use local "bindingEnergy()" to call through.
// 20100623  M. Kelsey -- Undo G4NucleiModel change from 0617.  Does not work
//		properly across multiple interactions.
// 20100627  M. Kelsey -- Protect recoil nucleus energy from floating roundoff
//		by setting small +ve or -ve values to zero.
// 20100701  M. Kelsey -- Let excitation energy be handled by G4InuclNuclei,
//		allow for ground-state recoil (goodCase == true for Eex==0.)
// 20100702  M. Kelsey -- Negative energy recoil should be rejected
// 20100706  D. Wright -- Copy "abandoned" cparticles to output list, copy
//		mesonic "excitons" to output list; should be absorbed, fix up
//		diagnostic messages.
// 20100713  M. Kelsey -- Add more diagnostics for Dennis' changes.
// 20100714  M. Kelsey -- Switch to new G4CascadeColliderBase class, remove
//		sanity check on afin/zfin (not valid).
// 20100715  M. Kelsey -- Add diagnostic for ekin_in vs. actual ekin; reduce
//		KE across Coulomb barrier.  Rearrange end-of-loop if blocks,
//		add conservation check at end.
// 20100716  M. Kelsey -- Eliminate inter_case; use base-class functionality.
//		Add minimum-fragment requirement for recoil, in order to
//		allow for momentum balancing
// 20100720  M. Kelsey -- Make EPCollider pointer member
// 20100721  M. Kelsey -- Turn on conservation checks unconditionally (override
//		new G4CASCADE_CHECK_ECONS setting
// 20100722  M. Kelsey -- Move cascade output buffers to .hh file
// 20100728  M. Kelsey -- Make G4NucleiModel data member for persistence,
//		delete colliders in destructor
// 20100906  M. Kelsey -- Hide "non-physical fragment" behind verbose flag
// 20100907  M. Kelsey -- Add makeResidualFragment function to create object
// 20100909  M. Kelsey -- Remove all local "fragment" stuff, use RecoilMaker.
//		move goodCase() to RecoilMaker.
// 20100910  M. Kelsey -- Use RecoilMaker::makeRecoilFragment().
// 20100915  M. Kelsey -- Define functions to deal with trapped particles,
//		move the exciton container to a data member
// 20100916  M. Kelsey -- Put decay photons directly onto output list
// 20100921  M. Kelsey -- Migrate to RecoilMaker::makeRecoilNuclei().
// 20100924  M. Kelsey -- Minor shuffling of post-cascade recoil building.
//		Create G4Fragment for recoil and store in output.
// 20110131  M. Kelsey -- Move "momentum_in" calculation inside verbosity
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110224  M. Kelsey -- Add ::rescatter() function which takes a list of
//		pre-existing secondaries as input.  Split ::collide() into
//		separate utility functions.  Move cascade parameters to static
//		data members.  Add setVerboseLevel().
// 20110302  M. Kelsey -- Move G4NucleiModel::printModel() call to G4NucleiModel
// 20110303  M. Kelsey -- Add more cascade functions to support rescattering
// 20110304  M. Kelsey -- Get original Propagate() arguments here in rescatter()
//		and convert to particles, nuclei and G4NucleiModel state.
// 20110308  M. Kelsey -- Don't put recoiling fragment onto output list any more
// 20110308  M. Kelsey -- Decay unstable hadrons from pre-cascade, use daughters
// 20110324  M. Kelsey -- Get locations of hit nuclei in ::rescatter(), pass
//		to G4NucleiModel::reset().
// 20110404  M. Kelsey -- Reduce maximum number of retries to 100, reflection
//		cut to 50.
// 20110721  M. Kelsey -- Put unusable pre-cascade particles directly on output,
//		do not decay.
// 20110722  M. Kelsey -- Deprecate "output_particles" list in favor of using
//		output directly (will help with pre-cascade issues).
// 20110801  M. Kelsey -- Use G4Inucl(Particle)::fill() functions to reduce
//		creation of temporaries.  Add local target buffer for
//		rescattering, to avoid memory leak.
// 20110808  M. Kelsey -- Pass buffer to generateParticleFate() to avoid copy
// 20110919  M. Kelsey -- Add optional final-state clustering, controlled (for
//		now) with compiler flag G4CASCADE_DO_COALESCENCE
// 20110922  M. Kelsey -- Follow migrations G4InuclParticle::print(ostream&)
//		and G4CascadParticle::print(ostream&); drop Q,B printing
// 20110926  M. Kelsey -- Replace compiler flag with one-time envvar in ctor
//		for final-state clustering.
// 20111003  M. Kelsey -- Prepare for gamma-N interactions by checking for
//		final-state tables instead of particle "isPhoton()"
// 20120521  A. Ribon -- Specify mass when decay trapped particle.
// 20120822  M. Kelsey -- Move envvars to G4CascadeParameters.
// 20121205  M. Kelsey -- In processSecondary(), set generation to 1, as these
//		particles are not true projectiles, but already embedded.
// 20130304  M. Kelsey -- Use new G4CascadeHistory to dump cascade structure
// 20140310  M. Kelsey -- (Bug #1584) Release memory allocated by DecayIt()
// 20140409  M. Kelsey -- Use const G4ParticleDefinition* everywhere
// 20141204  M. Kelsey -- Add function to test for non-interacting particles,
//		move those directly to output without propagating
// 20150608  M. Kelsey -- Label all while loops as terminating.
// 20150619  M. Kelsey -- Replace std::exp with G4Exp

#include <algorithm>

#include "G4IntraNucleiCascader.hh"
#include "G4SystemOfUnits.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeCoalescence.hh"
#include "G4CascadeHistory.hh"
#include "G4CascadeParameters.hh"
#include "G4CascadeRecoilMaker.hh"
#include "G4CascadParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4DecayProducts.hh"
#include "G4DecayTable.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4ExitonConfiguration.hh"
#include "G4Exp.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticleNames.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4LorentzConvertor.hh"
#include "G4Neutron.hh"
#include "G4NucleiModel.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4Proton.hh"
#include "G4V3DNucleus.hh"
#include "Randomize.hh"

using namespace G4InuclParticleNames;
using namespace G4InuclSpecialFunctions;


// Configuration parameters for cascade production
const G4int    G4IntraNucleiCascader::itry_max = 100;
const G4int    G4IntraNucleiCascader::reflection_cut = 50;
const G4double G4IntraNucleiCascader::small_ekin = 0.001*MeV;
const G4double G4IntraNucleiCascader::quasielast_cut = 1*MeV;


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4IntraNucleiCascader::G4IntraNucleiCascader()
  : G4CascadeColliderBase("G4IntraNucleiCascader"), model(new G4NucleiModel),
    theElementaryParticleCollider(new G4ElementaryParticleCollider),
    theRecoilMaker(new G4CascadeRecoilMaker), theClusterMaker(0),
    theCascadeHistory(0), tnuclei(0), bnuclei(0), bparticle(0),
    minimum_recoil_A(0.), coulombBarrier(0.),
    nucleusTarget(new G4InuclNuclei),
    protonTarget(new G4InuclElementaryParticle) {
  if (G4CascadeParameters::doCoalescence())
    theClusterMaker = new G4CascadeCoalescence;

  if (G4CascadeParameters::showHistory())
    theCascadeHistory = new G4CascadeHistory;
}

G4IntraNucleiCascader::~G4IntraNucleiCascader() {
  delete model;
  delete theElementaryParticleCollider;
  delete theRecoilMaker;
  delete theClusterMaker;
  delete theCascadeHistory;
  delete nucleusTarget;
  delete protonTarget;
}

void G4IntraNucleiCascader::setVerboseLevel(G4int verbose) {
  G4CascadeColliderBase::setVerboseLevel(verbose);
  model->setVerboseLevel(verbose);
  theElementaryParticleCollider->setVerboseLevel(verbose);
  theRecoilMaker->setVerboseLevel(verbose);

  // Optional functionality
  if (theClusterMaker) theClusterMaker->setVerboseLevel(verbose);
  if (theCascadeHistory) theCascadeHistory->setVerboseLevel(verbose);
}



void G4IntraNucleiCascader::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& globalOutput) {
  if (verboseLevel) G4cout << " >>> G4IntraNucleiCascader::collide " << G4endl;

  if (!initialize(bullet, target)) return;	// Load buffers and drivers

  G4int itry = 0;
  do {				/* Loop checking 08.06.2015 MHK */
    newCascade(++itry);
    setupCascade();
    generateCascade();
  } while (!finishCascade() && itry<itry_max);

  // Report full structure of final cascade if requested
  if (theCascadeHistory) theCascadeHistory->Print(G4cout);

  finalize(itry, bullet, target, globalOutput);
}

// For use with Propagate to preload a set of secondaries
// FIXME:  So far, we don't have any bullet information from Propagate!

void G4IntraNucleiCascader::rescatter(G4InuclParticle* bullet,
				      G4KineticTrackVector* theSecondaries,
				      G4V3DNucleus* theNucleus,
				      G4CollisionOutput& globalOutput) {
  if (verboseLevel)
    G4cout << " >>> G4IntraNucleiCascader::rescatter " << G4endl;

  G4InuclParticle* target = createTarget(theNucleus);
  if (!initialize(bullet, target)) return;	// Load buffers and drivers

  G4int itry = 0;
  do {				/* Loop checking 08.06.2015 MHK */
    newCascade(++itry);
    preloadCascade(theNucleus, theSecondaries);
    generateCascade();
  } while (!finishCascade() && itry<itry_max);

  // Report full structure of final cascade if requested
  if (theCascadeHistory) theCascadeHistory->Print(G4cout);

  finalize(itry, bullet, target, globalOutput);
}


G4bool G4IntraNucleiCascader::initialize(G4InuclParticle* bullet,
					 G4InuclParticle* target) {
  if (verboseLevel>1)
    G4cout << " >>> G4IntraNucleiCascader::initialize " << G4endl;
  
  // Configure processing modules
  theRecoilMaker->setTolerance(small_ekin);

  interCase.set(bullet,target);		// Classify collision type

  if (verboseLevel > 3) {
    G4cout << *interCase.getBullet() << G4endl
	   << *interCase.getTarget() << G4endl;
  }
  
  // Bullet may be nucleus or simple particle
  bnuclei = dynamic_cast<G4InuclNuclei*>(interCase.getBullet());
  bparticle = dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());
  
  if (!bnuclei && !bparticle) {
    G4cerr << " G4IntraNucleiCascader: projectile is not a valid particle."
	   << G4endl;
    return false;
  }
  
  // Target _must_ be nucleus
  tnuclei = dynamic_cast<G4InuclNuclei*>(interCase.getTarget());
  if (!tnuclei) {
    if (verboseLevel)
      G4cerr << " Target is not a nucleus.  Abandoning." << G4endl;
    return false;
  }
  
  model->generateModel(tnuclei);
  coulombBarrier = 0.00126*tnuclei->getZ() / (1.+G4cbrt(tnuclei->getA()));

  // Energy/momentum conservation usually requires a recoiling nuclear fragment
  // This cut will be increased on each "itry" if momentum could not balance.
  minimum_recoil_A = 0.;
  
  if (verboseLevel > 3) {
    G4LorentzVector momentum_in = bullet->getMomentum() + target->getMomentum();
    G4cout << " intitial momentum  E " << momentum_in.e() << " Px "
	   << momentum_in.x() << " Py " << momentum_in.y() << " Pz "
	   << momentum_in.z() << G4endl;
  }
  
  return true;
}

// Re-initialize buffers for new attempt at cascade

void G4IntraNucleiCascader::newCascade(G4int itry) {
  if (verboseLevel > 1) {
    G4cout << " IntraNucleiCascader itry " << itry << " inter_case "
	   << interCase.code() << G4endl;
  }
  
  model->reset();    			// Start new cascade process
  output.reset();
  new_cascad_particles.clear();
  theExitonConfiguration.clear();
  cascad_particles.clear();		// List of initial secondaries

  if (theCascadeHistory) theCascadeHistory->Clear();
}


// Load initial cascade using nuclear-model calculations

void G4IntraNucleiCascader::setupCascade() {
  if (verboseLevel > 1)
    G4cout << " >>> G4IntraNucleiCascader::setupCascade" << G4endl;

  if (interCase.hadNucleus()) { 		// particle with nuclei
    if (verboseLevel > 3)
      G4cout << " bparticle charge " << bparticle->getCharge()
	     << " baryon number " << bparticle->baryon() << G4endl;
    
    cascad_particles.push_back(model->initializeCascad(bparticle));
  } else {				// nuclei with nuclei
    G4int ab = bnuclei->getA();
    G4int zb = bnuclei->getZ();
    
    G4NucleiModel::modelLists all_particles;    // Buffer to receive lists
    model->initializeCascad(bnuclei, tnuclei, all_particles);
    
    cascad_particles = all_particles.first;
    output.addOutgoingParticles(all_particles.second);
    
    if (cascad_particles.size() == 0) { // compound nuclei
      G4int i;
      
      for (i = 0; i < ab; i++) {
	G4int knd = i < zb ? 1 : 2;
	theExitonConfiguration.incrementQP(knd);
      };
      
      G4int ihn = G4int(2 * (ab-zb) * inuclRndm() + 0.5);
      G4int ihz = G4int(2 * zb * inuclRndm() + 0.5);
      
      for (i = 0; i < ihn; i++) theExitonConfiguration.incrementHoles(2);
      for (i = 0; i < ihz; i++) theExitonConfiguration.incrementHoles(1);
    }
  }	// if (interCase ...
}


// Generate one possible cascade (all secondaries, etc.)

void G4IntraNucleiCascader::generateCascade() {
  if (verboseLevel>1) G4cout << " generateCascade " << G4endl;

  /* Loop checking 08.06.2015 MHK */
  G4int iloop = 0;
  while (!cascad_particles.empty() && !model->empty()) {
    iloop++;
    
    if (verboseLevel > 2) {
      G4cout << " Iteration " << iloop << ": Number of cparticles "
	     << cascad_particles.size() << " last one: \n"
	     << cascad_particles.back() << G4endl;
    }

    // Record incident particle first, to get history ID
    if (theCascadeHistory) {
      theCascadeHistory->AddEntry(cascad_particles.back());
      if (verboseLevel > 2) {
	G4cout << " active cparticle got history ID "
	       << cascad_particles.back().getHistoryId() << G4endl;
      }
    }

    // If non-interacting particle, move directly to output
    if (!particleCanInteract(cascad_particles.back())) {
      if (verboseLevel > 2)
	G4cout << " particle is non-interacting; moving to output" << G4endl;

      output.addOutgoingParticle(cascad_particles.back().getParticle());
      cascad_particles.pop_back();
      continue;
    }

    // Generate interaction with nucleon
    model->generateParticleFate(cascad_particles.back(),
				theElementaryParticleCollider,
				new_cascad_particles);

    // Record interaction for later reporting (if desired)
    if (theCascadeHistory && new_cascad_particles.size()>1) 
      theCascadeHistory->AddVertex(cascad_particles.back(), new_cascad_particles);

    if (verboseLevel > 2) {
      G4cout << " After generate fate: New particles "
	     << new_cascad_particles.size() << G4endl
	     << " Discarding last cparticle from list " << G4endl;
    }
    
    cascad_particles.pop_back();
    
    // handle the result of a new step
    
    if (new_cascad_particles.size() == 1) { // last particle goes without interaction
      const G4CascadParticle& currentCParticle = new_cascad_particles[0];
      
      if (model->stillInside(currentCParticle)) {
	if (verboseLevel > 3)
	  G4cout << " particle still inside nucleus " << G4endl;
	
	if (currentCParticle.getNumberOfReflections() < reflection_cut &&
	    model->worthToPropagate(currentCParticle)) {
	  if (verboseLevel > 3) G4cout << " continue reflections " << G4endl;
	  cascad_particles.push_back(currentCParticle);
	} else {
	  processTrappedParticle(currentCParticle);
	}	// reflection or exciton
	
      } else { // particle about to leave nucleus - check for Coulomb barrier
	if (verboseLevel > 3) G4cout << " possible escape " << G4endl;
	
	const G4InuclElementaryParticle& currentParticle =
	  currentCParticle.getParticle();
	
	G4double KE = currentParticle.getKineticEnergy();
	G4double mass = currentParticle.getMass();
	G4double Q = currentParticle.getCharge();
	
	if (verboseLevel > 3)
	  G4cout << " KE " << KE << " barrier " << Q*coulombBarrier << G4endl;
	
	if (KE < Q*coulombBarrier) {
	  // Calculate barrier penetration
	  G4double CBP = 0.0; 
	  
	  // if (KE > 0.0001) CBP = std::exp(-0.00126*tnuclei->getZ()*0.25*
	  //   (1./KE - 1./coulombBarrier));
	  if (KE > 0.0001) CBP = G4Exp(-0.0181*0.5*tnuclei->getZ()*
				       (1./KE - 1./coulombBarrier)*
				       std::sqrt(mass*(coulombBarrier-KE)) );
	  
	  if (G4UniformRand() < CBP) {
	    if (verboseLevel > 3) 
	      G4cout << " tunneled\n" << currentParticle << G4endl;

	    // Tunnelling through barrier leaves KE unchanged
	    output.addOutgoingParticle(currentParticle);
	  } else {
	    processTrappedParticle(currentCParticle);
	  }
	} else {
	  output.addOutgoingParticle(currentParticle);
	  
	  if (verboseLevel > 3)
	    G4cout << " Goes out\n" << output.getOutgoingParticles().back()
		   << G4endl;
	}
      } 
    } else { // interaction 
      if (verboseLevel > 3)
	G4cout << " interacted, adding new to list " << G4endl;
      
      cascad_particles.insert(cascad_particles.end(),
			      new_cascad_particles.begin(),
			      new_cascad_particles.end());
      
      std::pair<G4int, G4int> holes = model->getTypesOfNucleonsInvolved();
      if (verboseLevel > 3)
	G4cout << " adding new exciton holes " << holes.first << ","
	       << holes.second << G4endl;
      
      theExitonConfiguration.incrementHoles(holes.first);
      
      if (holes.second > 0)
	theExitonConfiguration.incrementHoles(holes.second);
    }		// if (new_cascad_particles ...
    
    // Evaluate nuclear residue
    theRecoilMaker->collide(interCase.getBullet(), interCase.getTarget(),
			    output, cascad_particles);
    
    G4double aresid = theRecoilMaker->getRecoilA();
    if (verboseLevel > 2) {
      G4cout << " cparticles remaining " << cascad_particles.size()
	     << " nucleus (model) has "
	     << model->getNumberOfNeutrons() << " n, "
	     << model->getNumberOfProtons() << " p "
	     << " residual fragment A " << aresid << G4endl;
    }
    
    if (aresid <= minimum_recoil_A) return;	// Must have minimum fragment
  }	// while cascade-list and model
}


// Conslidate results of cascade and evaluate success

G4bool G4IntraNucleiCascader::finishCascade() {
  if (verboseLevel > 1)
    G4cout << " >>> G4IntraNucleiCascader::finishCascade ?" << G4endl;

  // Add left-over cascade particles to output
  output.addOutgoingParticles(cascad_particles);
  cascad_particles.clear();

  // Cascade is finished. Check if it's OK.
  if (verboseLevel>3) {
    G4cout << " G4IntraNucleiCascader finished" << G4endl;
    output.printCollisionOutput();
  }

  // Apply cluster coalesence model to produce light ions
  if (theClusterMaker) {
    theClusterMaker->setVerboseLevel(verboseLevel);
    theClusterMaker->FindClusters(output);

    // Update recoil fragment after generating light ions
    if (verboseLevel>3) G4cout << " Recomputing recoil fragment" << G4endl;
    theRecoilMaker->collide(interCase.getBullet(), interCase.getTarget(),
			    output);
    if (verboseLevel>3) {
      G4cout << " After cluster coalescence" << G4endl;
      output.printCollisionOutput();
    }
  }

  // Use last created recoil fragment instead of re-constructing
  G4int afin = theRecoilMaker->getRecoilA();
  G4int zfin = theRecoilMaker->getRecoilZ();

  // FIXME:  Should we deal with unbalanced (0,0) case before rejecting?

  // Sanity check before proceeding
  if (!theRecoilMaker->goodFragment() && !theRecoilMaker->wholeEvent()) {
    if (verboseLevel > 1)
      G4cerr << " Recoil nucleus is not physical: A=" << afin << " Z="
	     << zfin << G4endl;
    return false;				// Discard event and try again
  }
  
  const G4LorentzVector& presid = theRecoilMaker->getRecoilMomentum();
  
  if (verboseLevel > 1) {
    G4cout << "  afin " << afin << " zfin " << zfin <<  G4endl;
  }
  
  if (afin == 0) return true;		// Whole event fragmented, exit
  
  if (afin == 1) {			// Add bare nucleon to particle list
    G4int last_type = (zfin==1) ? 1 : 2;	// proton=1, neutron=2
    
    G4double mass = G4InuclElementaryParticle::getParticleMass(last_type);
    G4double mres = presid.m();
    
    // Check for sensible kinematics
    if (mres-mass < -small_ekin) {		// Insufficient recoil energy
      if (verboseLevel > 2) G4cerr << " unphysical recoil nucleon" << G4endl;
      return false;
    }
    
    if (mres-mass > small_ekin) {		// Too much extra energy
      if (verboseLevel > 2)
	G4cerr << " extra energy with recoil nucleon" << G4endl;
      
      // FIXME:  For now, we add the nucleon as unbalanced, and let
      //	   "SetOnShell" fudge things.  This should be abandoned.
    }
    
    G4InuclElementaryParticle last_particle(presid, last_type, 
					    G4InuclParticle::INCascader);
    
    if (verboseLevel > 3) {
      G4cout << " adding recoiling nucleon to output list\n"
	     << last_particle  << G4endl;
    }
    
    output.addOutgoingParticle(last_particle);

    // Update recoil to include residual nucleon
    theRecoilMaker->collide(interCase.getBullet(), interCase.getTarget(),
			    output);
  }
  
  // Process recoil fragment for consistency, exit or reject
  if (output.numberOfOutgoingParticles() == 1) {
    G4double Eex = theRecoilMaker->getRecoilExcitation();
    if (std::abs(Eex) < quasielast_cut) {
      if (verboseLevel > 3) {
	G4cout << " quasi-elastic scatter with " << Eex << " MeV recoil"
	       << G4endl;
      }
      
      theRecoilMaker->setRecoilExcitation(Eex=0.);
      if (verboseLevel > 3) {
	G4cout << " Eex reset to " << theRecoilMaker->getRecoilExcitation()
	       << G4endl;
      }
    }
  }
  
  if (theRecoilMaker->goodNucleus()) {
    theRecoilMaker->addExcitonConfiguration(theExitonConfiguration);
    
    G4Fragment* recoilFrag = theRecoilMaker->makeRecoilFragment();
    if (!recoilFrag) {
      G4cerr << "Got null pointer for recoil fragment!" << G4endl;
      return false;
    }
    
    if (verboseLevel > 2)
      G4cout << " adding recoil fragment to output list" << G4endl;

    output.addRecoilFragment(*recoilFrag);
  }
  
  // Put final-state particles in "leading order" for return
  std::vector<G4InuclElementaryParticle>& opart = output.getOutgoingParticles();
  std::sort(opart.begin(), opart.end(), G4ParticleLargerEkin());
  
  // Adjust final state to balance momentum and energy if necessary
  if (theRecoilMaker->wholeEvent() || theRecoilMaker->goodNucleus()) {
    output.setVerboseLevel(verboseLevel);
    output.setOnShell(interCase.getBullet(), interCase.getTarget());
    output.setVerboseLevel(0);
    
    if (output.acceptable()) return true;
    else if (verboseLevel>2) G4cerr << " Cascade setOnShell failed." << G4endl;
  }

  // Cascade not physically reasonable
  if (afin <= minimum_recoil_A && minimum_recoil_A < tnuclei->getA()) {
    ++minimum_recoil_A;
    if (verboseLevel > 3) {
      G4cout << " minimum recoil fragment increased to A " << minimum_recoil_A
	     << G4endl;
    }
  }

  if (verboseLevel>2) G4cerr << " Cascade failed.  Retrying..." << G4endl;
  return false;
}


// Transfer finished cascade to return buffer

void 
G4IntraNucleiCascader::finalize(G4int itry, G4InuclParticle* bullet,
				G4InuclParticle* target,
				G4CollisionOutput& globalOutput) {
  if (itry >= itry_max) {
    if (verboseLevel) {
      G4cout << " IntraNucleiCascader-> no inelastic interaction after "
	     << itry << " attempts " << G4endl;
    }

    output.trivialise(bullet, target);
  } else if (verboseLevel) {
    G4cout << " IntraNucleiCascader output after trials " << itry << G4endl;
  }
  
  // Copy final generated cascade to output buffer for return
  globalOutput.add(output);
}


// Create simple nucleus from rescattering target

G4InuclParticle* 
G4IntraNucleiCascader::createTarget(G4V3DNucleus* theNucleus) {
  G4int theNucleusA = theNucleus->GetMassNumber();
  G4int theNucleusZ = theNucleus->GetCharge();
  
  if (theNucleusA > 1) {
    if (!nucleusTarget) nucleusTarget = new G4InuclNuclei;	// Just in case
    nucleusTarget->fill(0., theNucleusA, theNucleusZ, 0.);
    return nucleusTarget;
  } else {
    if (!protonTarget) protonTarget = new G4InuclElementaryParticle;
    protonTarget->fill(0., (theNucleusZ==1)?proton:neutron);
    return protonTarget;
  }

  return 0;		// Can never actually get here
}

// Copy existing (rescattering) cascade for propagation

void 
G4IntraNucleiCascader::preloadCascade(G4V3DNucleus* theNucleus,
				      G4KineticTrackVector* theSecondaries) {
  if (verboseLevel > 1)
    G4cout << " >>> G4IntraNucleiCascader::preloadCascade" << G4endl;

  copyWoundedNucleus(theNucleus);	// Update interacted nucleon counts
  copySecondaries(theSecondaries);	// Copy original to internal list
}

void G4IntraNucleiCascader::copyWoundedNucleus(G4V3DNucleus* theNucleus) {
  if (verboseLevel > 1)
    G4cout << " >>> G4IntraNucleiCascader::copyWoundedNucleus" << G4endl;

  // Loop over nucleons and count hits as exciton holes
  theExitonConfiguration.clear();
  hitNucleons.clear();
  if (theNucleus->StartLoop()) {
    G4Nucleon* nucl = 0;
    G4int nuclType = 0;
    /* Loop checking 08.06.2015 MHK */
    while ((nucl = theNucleus->GetNextNucleon())) {
      if (nucl->AreYouHit()) {	// Found previously interacted nucleon
	nuclType = G4InuclElementaryParticle::type(nucl->GetParticleType());
	theExitonConfiguration.incrementHoles(nuclType);
	hitNucleons.push_back(nucl->GetPosition());
      }
    }
  }

  if (verboseLevel > 3)
    G4cout << " nucleus has " << theExitonConfiguration.neutronHoles
	   << " neutrons hit, " << theExitonConfiguration.protonHoles
	   << " protons hit" << G4endl;

  // Preload nuclear model with confirmed hits, including locations
  model->reset(theExitonConfiguration.neutronHoles,
	       theExitonConfiguration.protonHoles, &hitNucleons);
}

void 
G4IntraNucleiCascader::copySecondaries(G4KineticTrackVector* secondaries) {
  if (verboseLevel > 1)
    G4cout << " >>> G4IntraNucleiCascader::copySecondaries" << G4endl;

  for (size_t i=0; i<secondaries->size(); i++) {
    if (verboseLevel > 3) G4cout << " processing secondary " << i << G4endl;

    processSecondary((*secondaries)[i]);      	// Copy to cascade or to output
  }

  // Sort list of secondaries to put leading particle first
  std::sort(cascad_particles.begin(), cascad_particles.end(),
	    G4ParticleLargerEkin());

  if (verboseLevel > 2) {
    G4cout << " Original list of " << secondaries->size() << " secondaries"
	   << " produced " << cascad_particles.size() << " cascade, "
	   << output.numberOfOutgoingParticles() << " released particles, "
	   << output.numberOfOutgoingNuclei() << " fragments" << G4endl;
  }
}


// Convert from pre-cascade secondary to local version

void G4IntraNucleiCascader::processSecondary(const G4KineticTrack* ktrack) {
  if (!ktrack) return;			// Sanity check

  // Get particle type to determine whether to keep or release
  const G4ParticleDefinition* kpd = ktrack->GetDefinition();
  if (!kpd) return;

  G4int ktype = G4InuclElementaryParticle::type(kpd);
  if (!ktype) {
    releaseSecondary(ktrack);
    return;
  }

  if (verboseLevel > 1) {
    G4cout << " >>> G4IntraNucleiCascader::processSecondary "
	   << kpd->GetParticleName() << G4endl;
  }

  // Allocate next local particle in buffer and fill
  cascad_particles.resize(cascad_particles.size()+1);	// Like push_back();
  G4CascadParticle& cpart = cascad_particles.back();

  // Convert momentum to Bertini internal units
  cpart.getParticle().fill(ktrack->Get4Momentum()/GeV, ktype);
  cpart.setGeneration(1);
  cpart.setMovingInsideNuclei();
  cpart.initializePath(0);

  // Convert position units to Bertini's internal scale
  G4ThreeVector cpos = ktrack->GetPosition()/model->getRadiusUnits();

  cpart.updatePosition(cpos);
  cpart.updateZone(model->getZone(cpos.mag()));

  if (verboseLevel > 2)
    G4cout << " Created cascade particle \n" << cpart << G4endl;
}


// Transfer unusable pre-cascade secondaries directly to output

void G4IntraNucleiCascader::releaseSecondary(const G4KineticTrack* ktrack) {
  const G4ParticleDefinition* kpd = ktrack->GetDefinition();

  if (verboseLevel > 1) {
    G4cout << " >>> G4IntraNucleiCascader::releaseSecondary "
	   << kpd->GetParticleName() << G4endl;
  }

  // Convert light ion into nucleus on fragment list
  if (dynamic_cast<const G4Ions*>(kpd)) {
    // Use resize() and fill() to avoid memory churn
    output.getOutgoingNuclei().resize(output.numberOfOutgoingNuclei()+1);
    G4InuclNuclei& inucl = output.getOutgoingNuclei().back();

    inucl.fill(ktrack->Get4Momentum()/GeV,
	       kpd->GetAtomicMass(), kpd->GetAtomicNumber());
    if (verboseLevel > 2)
      G4cout << " Created pre-cascade fragment\n" << inucl << G4endl;
  } else {
    // Use resize() and fill() to avoid memory churn
    output.getOutgoingParticles().resize(output.numberOfOutgoingParticles()+1);
    G4InuclElementaryParticle& ipart = output.getOutgoingParticles().back();

    // SPECIAL:  Use G4PartDef directly, allowing unknown type code
    ipart.fill(ktrack->Get4Momentum()/GeV, ktrack->GetDefinition());
    if (verboseLevel > 2)
      G4cout << " Created invalid pre-cascade particle\n" << ipart << G4endl;
  }
}

  
// Convert particles which cannot escape into excitons (or eject/decay them)

void G4IntraNucleiCascader::
processTrappedParticle(const G4CascadParticle& trapped) {
  const G4InuclElementaryParticle& trappedP = trapped.getParticle();

  G4int xtype = trappedP.type();
  if (verboseLevel > 3) G4cout << " exciton of type " << xtype << G4endl;
  
  if (trappedP.nucleon()) {	// normal exciton (proton or neutron)
    theExitonConfiguration.incrementQP(xtype);
    if (theCascadeHistory) theCascadeHistory->DropEntry(trapped);
    return;
  }

  if (trappedP.hyperon()) {	// Not nucleon, so must be hyperon
    decayTrappedParticle(trapped);
    if (theCascadeHistory) theCascadeHistory->DropEntry(trapped);
    return;
  }

  // non-standard exciton; release it
  // FIXME: this is a meson, so need to absorb it
  if (verboseLevel > 3) {
    G4cout << " non-standard should be absorbed, now released\n"
	   << trapped << G4endl;
  }
  
  output.addOutgoingParticle(trappedP);
}


// Decay unstable trapped particles, and add secondaries to processing list

void G4IntraNucleiCascader::
decayTrappedParticle(const G4CascadParticle& trapped) {
  if (verboseLevel > 3) 
    G4cout << " unstable must be decayed in flight" << G4endl;

  const G4InuclElementaryParticle& trappedP = trapped.getParticle();

  G4DecayTable* unstable = trappedP.getDefinition()->GetDecayTable();
  if (!unstable) {			// No decay table; cannot decay!
    if (verboseLevel > 3)
      G4cerr << " no decay table!  Releasing trapped particle" << G4endl;

    output.addOutgoingParticle(trappedP);
    return;
  }

  // Get secondaries from decay in particle's rest frame
  G4DecayProducts* daughters = unstable->SelectADecayChannel()->DecayIt( trappedP.getDefinition()->GetPDGMass() );
  if (!daughters) {			// No final state; cannot decay!
    if (verboseLevel > 3)
      G4cerr << " no daughters!  Releasing trapped particle" << G4endl;

    output.addOutgoingParticle(trappedP);
    return;
  }

  if (verboseLevel > 3)
    G4cout << " " << daughters->entries() << " decay daughters" << G4endl;

  // Convert secondaries to lab frame
  G4double decayEnergy = trappedP.getEnergy();
  G4ThreeVector decayDir = trappedP.getMomentum().vect().unit();
  daughters->Boost(decayEnergy, decayDir);

  // Put all the secondaries onto the list for propagation
  const G4ThreeVector& decayPos = trapped.getPosition();
  G4int zone = trapped.getCurrentZone();
  G4int gen = trapped.getGeneration()+1;

  for (G4int i=0; i<daughters->entries(); i++) {
    G4DynamicParticle* idaug = (*daughters)[i];

    G4InuclElementaryParticle idaugEP(*idaug, G4InuclParticle::INCascader);

    // Propagate hadronic secondaries with known interactions (tables)
    if (G4CascadeChannelTables::GetTable(idaugEP.type())) {
      if (verboseLevel > 3) G4cout << " propagating " << idaugEP << G4endl;
      cascad_particles.push_back(G4CascadParticle(idaugEP,decayPos,zone,0.,gen));
    } else {
      if (verboseLevel > 3) G4cout << " releasing " << idaugEP << G4endl;
      output.addOutgoingParticle(idaugEP);
    }
  }

  delete daughters;		// Clean up memory created by DecayIt()
}


// Test if particle is able to interact in nucleus

G4bool G4IntraNucleiCascader::
particleCanInteract(const G4CascadParticle& cpart) const {
  // If we have a lookup table for particle type on proton, it interacts
  return (0 != G4CascadeChannelTables::GetTable(cpart.getParticle().type()));
}
