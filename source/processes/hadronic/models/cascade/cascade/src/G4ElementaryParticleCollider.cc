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
// $Id$
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100316  D. Wright (restored by M. Kelsey) -- Replace original (incorrect)
//		pp, pn, nn 2-body to 2-body scattering angular distributions
//		with a new parametrization of elastic scattering data using
//		the sum of two exponentials.
// 20100319  M. Kelsey -- Use new generateWithRandomAngles for theta,phi stuff;
//		eliminate some unnecessary std::pow()
// 20100407  M. Kelsey -- Replace std::vector<>::resize(0) with ::clear()
//		Eliminate return-by-value std::vector<> by creating buffers
//		buffers for particles, momenta, and particle types.
//		The following functions now return void and are non-const:
//		  ::generateSCMfinalState()
//		  ::generateMomModules()
//		  ::generateStrangeChannelPartTypes()
//		  ::generateSCMpionAbsorption()
// 20100408  M. Kelsey -- Follow changes to G4*Sampler to pass particle_kinds
//		as input buffer.
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100428  M. Kelsey -- Use G4InuclParticleNames enum
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()"
// 20100507  M. Kelsey -- Rationalize multiplicity returns to be actual value
// 20100511  M. Kelsey -- Replace G4PionSampler and G4NucleonSampler with new
//		pi-N and N-N classes, reorganize if-cascades 
// 20100511  M. Kelsey -- Eliminate three residual random-angles blocks.
// 20100511  M. Kelsey -- Bug fix: pi-N two-body final states not correctly
//		tested for charge-exchange case.
// 20100512  M. Kelsey -- Rationalize multiplicity returns to be actual value
// 20100512  M. Kelsey -- Add some additional debugging messages for 2-to-2
// 20100512  M. Kelsey -- Replace "if (is==)" cascades with switch blocks.
//		Use G4CascadeInterpolator for angular distributions.
// 20100517  M. Kelsey -- Inherit from common base class, make arrays static
// 20100519  M. Kelsey -- Use G4InteractionCase to compute "is" values.
// 20100625  M. Kelsey -- Two bugs in n-body momentum, last particle recoil
// 20100713  M. Kelsey -- Bump collide start message up to verbose > 1
// 20100714  M. Kelsey -- Move conservation checking to base class
// 20100714  M. Kelsey -- Add sanity check for two-body final state, to ensure
//		that final state total mass is below etot_scm; also compute
//		kinematics without "rescaling" (which led to non-conservation)
// 20100726  M. Kelsey -- Move remaining std::vector<> buffers to .hh file
// 20100804  M. Kelsey -- Add printing of final-state tables, protected by
//		G4CASCADE_DEBUG_SAMPLER preprocessor flag
// 20101019  M. Kelsey -- CoVerity report: check dynamic_cast<> for null
// 20110214  M. Kelsey -- Follow G4InuclParticle::Model enumerator migration
// 20110718  V. Uzhinskiy -- Drop IL,IM reset for multiplicity-three collisions
// 20110718  M. Kelsey -- Use enum names in switch blocks (c.f. G4NucleiModel)
// 20110720  M. Kelsey -- Follow interface change for cross-section tables,
//		eliminating switch blocks.
// 20110806  M. Kelsey -- Pre-allocate buffers to avoid memory churn
// 20110922  M. Kelsey -- Follow G4InuclParticle::print(ostream&) migration
// 20110926  M. Kelsey -- Protect sampleCMcosFor2to2 from unphysical arguments
// 20111003  M. Kelsey -- Prepare for gamma-N interactions by checking for
//		final-state tables instead of particle "isPhoton()"
// 20111007  M. Kelsey -- Add gamma-N final-state tables to printFinalState
// 20111107  M. Kelsey -- In sampleCMmomentumFor2to2(), hide message about
//		unrecognized gamma-N initial state behind verbosity.
// 20111216  M. Kelsey -- Add diagnostics to generateMomModulesFor2toMany(),
//		defer allocation of buffer in generateSCMfinalState() so that
//		multiplicity failures return zero output, and can be trapped.
// 20120308  M. Kelsey -- Allow photons to interact with dibaryons (see
//		changes in G4NucleiModel).
// 20120608  M. Kelsey -- Fix variable-name "shadowing" compiler warnings.
// 20121206  M. Kelsey -- Add Omega to printFinalStateTables(), remove line
//		about "preliminary" gamma-N.
// 20130129  M. Kelsey -- Add static arrays and interpolators for two-body
//		angular distributions (addresses MT thread-local issue)
// 20130221  M. Kelsey -- Move two-body angular dist classes to factory
// 20130306  M. Kelsey -- Use particle-name enums in if-blocks; add comments
//		to sections of momentum-coefficient matrix; move final state
//		table printing to G4CascadeChannelTables.
// 20130307  M. Kelsey -- Reverse order of dimensions for rmn array
// 20130307  M. Kelsey -- Use new momentum generator factory instead of rmn
// 20130308  M. Kelsey -- Move 3-body angle calc to G4InuclSpecialFunctions.

#include "G4ElementaryParticleCollider.hh"
#include "G4CascadeChannel.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclParticleNames.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4MultiBodyMomentumDist.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4TwoBodyAngularDist.hh"
#include "G4VMultiBodyMomDst.hh"
#include "G4VTwoBodyAngDst.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cfloat>
#include <vector>

using namespace G4InuclParticleNames;
using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;


G4ElementaryParticleCollider::G4ElementaryParticleCollider()
  : G4CascadeColliderBase("G4ElementaryParticleCollider") {;}


void
G4ElementaryParticleCollider::collide(G4InuclParticle* bullet,
				      G4InuclParticle* target,
				      G4CollisionOutput& output) 
{
  if (verboseLevel > 1)
    G4cout << " >>> G4ElementaryParticleCollider::collide" << G4endl;

  if (!useEPCollider(bullet,target)) {		// Sanity check
    G4cerr << " ElementaryParticleCollider -> can collide only particle with particle " 
           << G4endl;
    return;
  }

#ifdef G4CASCADE_DEBUG_SAMPLER
  static G4bool doPrintTables = true;	// Once and only once per job
  if (doPrintTables) {
    G4CascadeChannelTables::Print();
    doPrintTables = false;
  }
#endif

  interCase.set(bullet, target);	// To identify kind of collision

  if (verboseLevel > 1) G4cout << *bullet << G4endl << *target << G4endl;

  G4InuclElementaryParticle* particle1 =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* particle2 =	
    dynamic_cast<G4InuclElementaryParticle*>(target);

  if (!particle1 || !particle2) {	// Redundant with useEPCollider()
    G4cerr << " ElementaryParticleCollider -> can only collide hadrons"
           << G4endl;
    return;
  }

  // Check for available interaction, or pion+dibaryon special case
  if (!G4CascadeChannelTables::GetTable(interCase.hadrons()) &&
      !particle1->quasi_deutron() && !particle2->quasi_deutron()) {
    G4cerr << " ElementaryParticleCollider -> cannot collide " 
	   << particle1->getDefinition()->GetParticleName() << " with "
           << particle2->getDefinition()->GetParticleName() << G4endl;
    return;
  }
  // Generate nucleon or pion collision with nucleon
  // or pion with quasi-deuteron

  if (particle1->nucleon() || particle2->nucleon()) { // ok
    G4LorentzConvertor convertToSCM;
    if(particle2->nucleon()) {
      convertToSCM.setBullet(particle1);
      convertToSCM.setTarget(particle2);
    } else {
      convertToSCM.setBullet(particle2);
      convertToSCM.setTarget(particle1);
    };

    convertToSCM.setVerbose(verboseLevel);

    convertToSCM.toTheCenterOfMass();
    G4double ekin = convertToSCM.getKinEnergyInTheTRS();
    G4double etot_scm = convertToSCM.getTotalSCMEnergy();
    G4double pscm = convertToSCM.getSCMMomentum();

    generateSCMfinalState(ekin, etot_scm, pscm, particle1, particle2,
			  &convertToSCM);

    if (particles.empty()) {	// No final state possible, pass bullet through
      if (verboseLevel) {
	G4cerr << " ElementaryParticleCollider -> failed to collide " 
	       << particle1->getMomModule() << " GeV/c " 
	       << particle1->getDefinition()->GetParticleName() << " with "
	       << particle2->getDefinition()->GetParticleName() << G4endl;
      }
    } else {			 // convert back to Lab
      G4LorentzVector mom;		// Buffer to avoid memory churn
      particleIterator ipart;
      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {	
	mom = convertToSCM.backToTheLab(ipart->getMomentum());
	ipart->setMomentum(mom); 
      };

      // Check conservation in multibody final state
      if (verboseLevel && !validateOutput(bullet, target, particles)) {
	G4cout << " incoming particles: \n" << *particle1 << G4endl
	       << *particle2 << G4endl
	       << " outgoing particles: " << G4endl;
	for(ipart = particles.begin(); ipart != particles.end(); ipart++)
	  G4cout << *ipart << G4endl;

	G4cout << " <<< Non-conservation in G4ElementaryParticleCollider"
	       << G4endl;
      }

      std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
      output.addOutgoingParticles(particles);
    }
  } else {	// neither particle is nucleon: pion on quasideuteron
    if (particle1->quasi_deutron() || particle2->quasi_deutron()) {
      if (particle1->pion() || particle2->pion() ||
	  particle1->isPhoton() || particle2->isPhoton()) {
	G4LorentzConvertor convertToSCM;
	if(particle2->quasi_deutron()) {	// Quasideuteron is target
	  convertToSCM.setBullet(particle1);
	  convertToSCM.setTarget(particle2);
	} else {
	  convertToSCM.setBullet(particle2);
	  convertToSCM.setTarget(particle1);
	}; 
	convertToSCM.toTheCenterOfMass(); 
	G4double etot_scm = convertToSCM.getTotalSCMEnergy();
	
	generateSCMpionAbsorption(etot_scm, particle1, particle2);

	if (particles.empty()) {	// Failed to generate final state
	  if (verboseLevel) {
	    G4cerr << " ElementaryParticleCollider -> failed to collide " 
		   << particle1->getMomModule() << " GeV/c " 
		   << particle1->getDefinition()->GetParticleName() << " with "
		   << particle2->getDefinition()->GetParticleName() << G4endl;
	  }
	} else {			// convert back to Lab
	  G4LorentzVector mom;	// Buffer to avoid memory churn
	  particleIterator ipart;
	  for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	    mom = convertToSCM.backToTheLab(ipart->getMomentum());
	    ipart->setMomentum(mom); 
	  };

	  validateOutput(bullet, target, particles);	// Check conservation

	  std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
	  output.addOutgoingParticles(particles);
	};
      } else {
	G4cerr << " ElementaryParticleCollider -> can only collide pions with dibaryons " 
	       << G4endl;
      };
    } else {
      G4cerr << " ElementaryParticleCollider -> can only collide something with nucleon or dibaryon " 
	     << G4endl;
    };
  };  
}


G4int 
G4ElementaryParticleCollider::generateMultiplicity(G4int is, 
						   G4double ekin) const 
{
  G4int mul = 0;

  const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(is);

  if (xsecTable) mul = xsecTable->getMultiplicity(ekin);
  else {
    G4cerr << " G4ElementaryParticleCollider: Unknown interaction channel "
	   << is << " - multiplicity not generated " << G4endl;
  }

  if(verboseLevel > 3){
    G4cout << " G4ElementaryParticleCollider::generateMultiplicity: "  
           << " multiplicity = " << mul << G4endl; 
  }

  return mul;
}

 
void
G4ElementaryParticleCollider::generateSCMfinalState(G4double ekin, 
		                     G4double etot_scm, 
		                     G4double pscm,
		                     G4InuclElementaryParticle* particle1,
		                     G4InuclElementaryParticle* particle2, 
	                             G4LorentzConvertor* toSCM) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMfinalState" 
           << G4endl;
  }

  const G4double ang_cut = 0.9999;
  const G4double difr_const = 0.3678794;   
  const G4int itry_max = 10;
  G4InuclElementaryParticle dummy;

  G4int type1 = particle1->type();
  G4int type2 = particle2->type();

  G4int is = type1 * type2;

  if(verboseLevel > 3){
    G4cout << " is " << is << G4endl;
  }

  G4int multiplicity = 0;
  G4bool generate = true;

  while (generate) {
    particles.clear();		// Initialize buffers for this event
    particle_kinds.clear();

    // Generate list of final-state particles
    multiplicity = generateMultiplicity(is, ekin);

    generateOutgoingPartTypes(is, multiplicity, ekin);
    if (particle_kinds.empty()) {
      if (verboseLevel > 3) {
	G4cout << " generateOutgoingPartTypes failed mult " << multiplicity
	       << G4endl;
      }
      continue;
    }

    if (multiplicity == 2) {
      // Identify charge or strangeness exchange (non-elastic scatter)
      G4int finaltype = particle_kinds[0]*particle_kinds[1];
      G4int kw = (finaltype != is) ? 2 : 1;

      G4double pmod = pscm;	// Elastic scattering preserves CM momentum

      if (kw == 2) {		// Non-elastic needs new CM momentum value
	G4double mone = dummy.getParticleMass(particle_kinds[0]);
	G4double mtwo = dummy.getParticleMass(particle_kinds[1]);

	if (etot_scm < mone+mtwo) {		// Can't produce final state
	  if (verboseLevel > 2) {
	    G4cerr << " bad final state " << particle_kinds[0]
		   << " , " << particle_kinds[1] << " etot_scm " << etot_scm
		   << " < mone+mtwo " << mone+mtwo << " , but ekin " << ekin
		   << G4endl;
	  }
	  continue;
	}

	G4double ecm_sq = etot_scm*etot_scm;
	G4double msumsq = mone+mtwo; msumsq *= msumsq;
	G4double mdifsq = mone-mtwo; mdifsq *= mdifsq;

	G4double a = (ecm_sq - msumsq) * (ecm_sq - mdifsq);

	pmod = std::sqrt(a)/(2.*etot_scm);
      }

      G4LorentzVector mom = sampleCMmomentumFor2to2(is, finaltype, kw, ekin, pmod);

      if (verboseLevel > 3) {
	G4cout << " Particle kinds = " << particle_kinds[0] << " , "
	       << particle_kinds[1] << G4endl
	       << " pscm " << pscm << " pmod " << pmod << G4endl
	       << " before rotation px " << mom.x() << " py " << mom.y()
	       << " pz " << mom.z() << G4endl;
      }

      mom = toSCM->rotate(mom); 

      if (verboseLevel > 3){
	G4cout << " after rotation px " << mom.x() << " py " << mom.y() <<
	  " pz " << mom.z() << G4endl;
      }
      G4LorentzVector mom1(-mom.vect(), mom.e());

      particles.resize(multiplicity);		// Preallocate buffer
      particles[0].fill(mom, particle_kinds[0], G4InuclParticle::EPCollider);
      particles[1].fill(mom1, particle_kinds[1], G4InuclParticle::EPCollider);
      generate = false;
    } else {			 // 2 -> many
      G4int itry = 0;
      G4bool bad = true;
      G4int knd_last = particle_kinds[multiplicity - 1];
      G4double mass_last = dummy.getParticleMass(knd_last);

      if (verboseLevel > 3){
	G4cout << " knd_last " << knd_last << " mass " << mass_last << G4endl;
      }

      while (bad && itry < itry_max) {
	itry++;

	if (verboseLevel > 3){
	  G4cout << " itry in generateSCMfinalState " << itry << G4endl;
	}

	generateMomModules(multiplicity, is, ekin, etot_scm);
	if (G4int(modules.size()) != multiplicity) {
	  if (verboseLevel > 3) {
	    G4cerr << " generateMomModule failed at mult " << multiplicity
		   << " ekin " << ekin << " etot_scm " << etot_scm << G4endl;
	  }
	  continue;
	}

	if (multiplicity == 3) { 
	  G4LorentzVector mom3 = 
	    particleSCMmomentumFor2to3(is, knd_last, ekin, modules[2]);
	  
	  mom3 = toSCM->rotate(mom3);
	  
	  // generate the momentum of first
	  G4double ct = -0.5 * (modules[2] * modules[2] + 
				modules[0] * modules[0] - 
				modules[1] * modules[1]) /
	    modules[2] / modules[0];   
	  
	  if(std::fabs(ct) < ang_cut) {
	    
	    if(verboseLevel > 2){
	      G4cout << " ok for mult " << multiplicity << G4endl;
	    }
	    
	    G4LorentzVector mom1 = generateWithFixedTheta(ct, modules[0]);
	    mom1 = toSCM->rotate(mom3, mom1);

	    // Second particle recoils off 1 & 3
	    G4LorentzVector mom2(etot_scm);
	    mom2 -= mom1+mom3;
	    
	    bad = false;
	    generate = false;

	    particles.resize(multiplicity);		// Preallocate buffer
   
	    particles[0].fill(mom1, particle_kinds[0], G4InuclParticle::EPCollider);
	    particles[1].fill(mom2, particle_kinds[1], G4InuclParticle::EPCollider);
	    particles[2].fill(mom3, particle_kinds[2], G4InuclParticle::EPCollider);
	  };
	} else { // multiplicity > 3
	  // generate first mult - 2 momentums
	  G4LorentzVector tot_mom;
	  scm_momentums.clear();
	  
	  for (G4int i = 0; i < multiplicity - 2; i++) {
	    G4double p0 = particle_kinds[i] < 3 ? 0.36 : 0.25;
	    G4double alf = 1.0 / p0 / (p0 - (modules[i] + p0) *
				       std::exp(-modules[i] / p0));
	    G4double st = 2.0;
	    G4int itry1 = 0;
	    
	    while (std::fabs(st) > ang_cut && itry1 < itry_max) {
	      itry1++;
	      G4double s1 = modules[i] * inuclRndm();
	      G4double s2 = alf * difr_const * p0 * inuclRndm();
	      
	      if(verboseLevel > 3){
		G4cout << " s1 * alf * std::exp(-s1 / p0) "
		       << s1 * alf * std::exp(-s1 / p0) 
		       << " s2 " << s2 << G4endl;
	      }
	      
	      if(s1 * alf * std::exp(-s1 / p0) > s2) st = s1 / modules[i];
	    };
	    
	    if(verboseLevel > 3){
	      G4cout << " itry1 " << itry1 << " i " << i << " st " << st
		     << G4endl;
	    }
	    
	    if(itry1 == itry_max) {
	      if(verboseLevel > 2){
		G4cout << " high energy angles generation: itry1 " << itry1
		       << G4endl;
	      }
	      
	      st = 0.5 * inuclRndm();
	    };

	    G4double ct = std::sqrt(1.0 - st * st);
	    if (inuclRndm() > 0.5) ct = -ct;
	    
	    G4LorentzVector mom = generateWithFixedTheta(ct,modules[i]);

	    tot_mom += mom;
	    
	    scm_momentums.push_back(mom);
	  }; 
	  
	  // handle last two
	  G4double tot_mod = tot_mom.rho(); 
	  G4double ct = -0.5 * (tot_mod * tot_mod + 
				modules[multiplicity - 2] * modules[multiplicity - 2] -
				modules[multiplicity - 1] * modules[multiplicity - 1]) / tot_mod /
	    modules[multiplicity - 2];  
	  
	  if (verboseLevel > 2){
	    G4cout << " ct last " << ct << G4endl;
	  }            
	  
	  if (std::fabs(ct) < ang_cut) {
	    G4int i(0);
	    for (i = 0; i < multiplicity - 2; i++) 
	      scm_momentums[i] = toSCM->rotate(scm_momentums[i]);
	    
	    tot_mom = toSCM->rotate(tot_mom);  
	    
	    G4LorentzVector mom = 
	      generateWithFixedTheta(ct, modules[multiplicity - 2]);
	    
	    mom = toSCM->rotate(tot_mom, mom);
	    scm_momentums.push_back(mom);

	    // Last particle recoils off everything else
	    G4LorentzVector mom1(etot_scm);
	    mom1 -= mom+tot_mom;
	    
	    scm_momentums.push_back(mom1);  
	    bad = false;
	    generate = false;
	    
	    if (verboseLevel > 2){
	      G4cout << " ok for mult " << multiplicity << G4endl;
	    }

	    particles.resize(multiplicity);		// Preallocate buffer
	    for (i = 0; i < multiplicity; i++) {
	      particles[i].fill(scm_momentums[i], particle_kinds[i],
				G4InuclParticle::EPCollider);
	    }
	  }	// |ct| < ang_cut
	}	// multiplicity > 3
      }		// while (bad && itry < itry_max)

      if (itry == itry_max) {
	if (verboseLevel > 2) {
	  G4cout << " cannot generate the distr. for mult " << multiplicity
		 << G4endl;
	}
	break;
      }
    }	// multiplicity > 2
  }	// while (generate) 

  if (verboseLevel > 3) {
    G4cout << " <<< G4ElementaryParticleCollider::generateSCMfinalState"
	   << G4endl;
  }

  return;	// Particles buffer filled
}

void 
G4ElementaryParticleCollider::generateMomModules(G4int mult, 
						 G4int is, 
						 G4double ekin, 
						 G4double etot_cm) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateMomModules" 
           << G4endl;
  }

  if (verboseLevel > 2){
    G4cout << " mult " << mult << " is " << is << " ekin " << ekin
	   << " etot_cm " << etot_cm << G4endl;
  }

  const G4int itry_max = 10;
  const G4double small = 1.e-10;
  G4InuclElementaryParticle dummy;
  G4int itry = 0;

  modules.clear();			// Initialize buffer for this attempt
  modules.resize(mult,0.);

  masses2.clear();
  masses2.resize(mult,0.);		// Allows direct [i] setting

  for (G4int i = 0; i < mult; i++) {
    G4double mass = dummy.getParticleMass(particle_kinds[i]);
    masses2[i] = mass * mass;
  };

  G4double mass_last = std::sqrt(masses2[mult - 1]);

  if (verboseLevel > 3){
    G4cout << " knd_last " << particle_kinds[mult - 1] << " mass_last " 
           << mass_last << G4endl;
  }

  while (itry < itry_max) {
    itry++;
    if(verboseLevel > 3){
      G4cout << " itry in generateMomModules " << itry << G4endl;
    }

    G4int ilast = -1;
    G4double eleft = etot_cm;

    for (G4int i = 0; i < mult - 1; i++) {
      G4double pmod = 
	getMomModuleFor2toMany(is, mult, particle_kinds[i], ekin);

      if (pmod < small) break;
      eleft -= std::sqrt(pmod * pmod + masses2[i]);

      if (verboseLevel > 3){
	G4cout << " kp " << particle_kinds[i] << " pmod " << pmod
	       << " mass2 " << masses2[i] << " eleft " << eleft
	       << "\n x1 " << eleft - mass_last << G4endl;
      }

      if (eleft <= mass_last) break;
      ilast++;
      modules[i] = pmod;
    };

    if (ilast == mult - 2) {
      G4double plast = eleft * eleft - masses2[mult - 1];
      if (verboseLevel > 2){
	G4cout << " plast ** 2 " << plast << G4endl;
      }

      if (plast > small) {
	plast = std::sqrt(plast);
	modules[mult - 1] = plast;      

	if (mult == 3) { 
	  if (satisfyTriangle(modules)) return;
	} else return;
      }
    }
  }	// while (itry < itry_max)

  if (verboseLevel > 2)
    G4cerr << " Unable to generate momenta for multiplicity " << mult << G4endl;

  modules.clear();		// Something went wrong, throw away partial
  return;    
}


G4bool 
G4ElementaryParticleCollider::satisfyTriangle(
			const std::vector<G4double>& pmod) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::satisfyTriangle" 
           << G4endl;
  }

  G4bool good = ( (pmod.size() != 3) ||
		  !(std::fabs(pmod[1] - pmod[2]) > pmod[0] || 
		    pmod[0] > pmod[1] + pmod[2] ||
		    std::fabs(pmod[0] - pmod[2]) > pmod[1] ||
		    pmod[1] > pmod[0] + pmod[2] ||
		    std::fabs(pmod[0] - pmod[1]) > pmod[2] ||
		    pmod[2] > pmod[1] + pmod[0]));

  return good;
}


void 
G4ElementaryParticleCollider::generateOutgoingPartTypes(G4int is, G4int mult,
							G4double ekin)
{
  particle_kinds.clear();	// Initialize buffer for generation

  const G4CascadeChannel* xsecTable = G4CascadeChannelTables::GetTable(is);

  if (xsecTable)
    xsecTable->getOutgoingParticleTypes(particle_kinds, mult, ekin);
  else {
    G4cerr << " G4ElementaryParticleCollider: Unknown interaction channel "
	   << is << " - outgoing kinds not generated " << G4endl;
  }

  return;
}


G4double 
G4ElementaryParticleCollider::getMomModuleFor2toMany(G4int is, G4int mult, 
					             G4int knd, 
					     	     G4double ekin) const 
{
  if (verboseLevel > 2) {
    G4cout << " >>> G4ElementaryParticleCollider::getMomModuleFor2toMany "
	   << " is " << is << " mult " << mult << " knd " << knd
	   << " ekin " << ekin << G4endl;
  }

  // Get generator for specified state
  const G4VMultiBodyMomDst* momDist = G4MultiBodyMomentumDist::GetDist(is,mult);
  if (!momDist && verboseLevel) {
    G4cerr << " G4ElementaryParticleCollider::getMomModuleFor2toMany:"
	   << " interaction is=" << is << " mult=" << mult << " not recognized "
	   << G4endl;
  }

  // Choose momentum modulus for state, or fraction of kinetic energy
  if (verboseLevel && momDist)
    const_cast<G4VMultiBodyMomDst*>(momDist)->setVerboseLevel(verboseLevel);

  G4double pmod = 0.;
  if (momDist) pmod = momDist->GetMomentum(knd, ekin);
  else pmod = inuclRndm() * ekin/mult;

  return pmod;
}


G4LorentzVector 
G4ElementaryParticleCollider::particleSCMmomentumFor2to3(
			   G4int is, 
			   G4int knd, 
			   G4double ekin, 
			   G4double pmod) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::particleSCMmomentumFor2to3" 
           << G4endl;
  }

  const G4int itry_max = 100;
  G4double ct = 2.0;
  G4int K = 2;
  G4int J = 1;

  if (is == pro*pro || is == pro*neu || is == neu*neu) K = 0;
  if (knd == proton || knd == neutron) J = 0;

  G4int itry = 0;

  while (std::fabs(ct) > 1.0 && itry < itry_max) {
    itry++;
    G4double Spow = randomInuclPowers(ekin, abn[K+J]);
    ct = 2.0*Spow - 1.0;
  }

  if (itry == itry_max) {	// No success, just throw flat distribution
    if (verboseLevel > 2) {
      G4cout << " particleSCMmomentumFor2to3 -> itry = itry_max " << itry
	     << G4endl;
    }

    ct = 2.0 * inuclRndm() - 1.0;
  };

  return generateWithFixedTheta(ct, pmod);
}


G4LorentzVector 
G4ElementaryParticleCollider::sampleCMmomentumFor2to2(G4int is, G4int fs, G4int kw, 
                                                      G4double ekin, 
			                              G4double pscm) const 
{
  if (verboseLevel > 3)
    G4cout << " >>> G4ElementaryParticleCollider::sampleCMmomentumFor2to2" 
	   << " is " << is << " kw " << kw << " ekin " << ekin << " pscm "
	   << pscm << G4endl;

  // Get generator for specified state
  const G4VTwoBodyAngDst* angDist = G4TwoBodyAngularDist::GetDist(is,fs,kw);
  if (!angDist && verboseLevel) {
    G4cerr << " G4ElementaryParticleCollider::sampleCMmomentumFor2to2:"
	   << " interaction is=" << is << " kw=" << kw << " not recognized "
	   << G4endl;
  }

  // Choose cos(theta) for state, or forward scatter if problem
  if (verboseLevel && angDist) 
    const_cast<G4VTwoBodyAngDst*>(angDist)->setVerboseLevel(verboseLevel);

  G4double ct = angDist ? angDist->GetCosTheta(ekin, pscm) : 1.;

  return generateWithFixedTheta(ct, pscm);
}


void
G4ElementaryParticleCollider::generateSCMpionAbsorption(G4double etot_scm,
			             G4InuclElementaryParticle* particle1,
			             G4InuclElementaryParticle* particle2) {
  if (verboseLevel > 3)
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMpionAbsorption" 
	   << G4endl;

  // generate nucleons momenta for pion or photon absorption
  // the nucleon distribution assumed to be isotropic in SCM

  particles.clear();		// Initialize buffers for this event
  particles.resize(2);

  particle_kinds.clear();

  G4int type1 = particle1->type();
  G4int type2 = particle2->type();

  // generate kinds

  if (type1 == pionPlus) {
    if (type2 == diproton) {		// pi+ + PP -> ? 
      G4cerr << " pion absorption: pi+ + PP -> ? " << G4endl;
      return;
    } else if (type2 == unboundPN) { 	// pi+ + PN -> PP
      particle_kinds.push_back(proton);
      particle_kinds.push_back(proton);
    } else if (type2 == dineutron) { 	// pi+ + NN -> PN
      particle_kinds.push_back(proton);
      particle_kinds.push_back(neutron);
    }
  } else if (type1 == pionMinus) { 
    if (type2 == diproton) {		// pi- + PP -> PN 
      particle_kinds.push_back(proton);
      particle_kinds.push_back(neutron);
    } else if (type2 == unboundPN) {	 // pi- + PN -> NN
      particle_kinds.push_back(neutron);
      particle_kinds.push_back(neutron);
    } else if (type2 == dineutron) {	// pi- + NN -> ?
      G4cerr << " pion absorption: pi- + NN -> ? " << G4endl;
      return;
    }
  } else if (type1 == pionZero || type1 == photon) {
    if (type2 == diproton) {		// pi0/gamma + PP -> PP 
      particle_kinds.push_back(proton);
      particle_kinds.push_back(proton);
    } else if (type2 == unboundPN) {	// pi0/gamma + PN -> PN
      particle_kinds.push_back(proton);
      particle_kinds.push_back(neutron);
    } else if (type2 == dineutron) {	// pi0/gamma + NN -> ?
      particle_kinds.push_back(neutron);
      particle_kinds.push_back(neutron);
    }
  }
    
  G4InuclElementaryParticle dummy;

  G4double mone = dummy.getParticleMass(particle_kinds[0]);
  G4double m1sq = mone*mone;

  G4double mtwo = dummy.getParticleMass(particle_kinds[1]);
  G4double m2sq = mtwo*mtwo;

  G4double a = 0.5 * (etot_scm * etot_scm - m1sq - m2sq);

  G4double pmod = std::sqrt((a * a - m1sq * m2sq) / (m1sq + m2sq + 2.0 * a));
  G4LorentzVector mom1 = generateWithRandomAngles(pmod, mone);
  G4LorentzVector mom2;
  mom2.setVectM(-mom1.vect(), mtwo);

  particles[0].fill(mom1, particle_kinds[0], G4InuclParticle::EPCollider);
  particles[1].fill(mom2, particle_kinds[1], G4InuclParticle::EPCollider);

  return;
}


// Parameter array for cos(theta) calculation of three body final states
// See particleSCMmomentumFor2to3() for variables mentioned below
//
// Outer:   [0] (K=0 NN initial state, J=0 outgoing N)
//          [1] (K=0 NN initial state, J=1 outgoing h,K,Y)
//          [2] (K=2 hN,KN,YN,gN initial state, J=0 outgoing N)
//          [3] (K=2 hN,KN,YN,gN initial state, J=1 outgoing h,K,Y)
//
// Middle:  blocks for powers of S^0..3
//
// Inner:   Powers of Ekin^0..3

const G4double G4ElementaryParticleCollider::abn[4][4][4] = {
  // -------- Initial state nucleon-nucleon, outgoing N --------
  { { 0.0856, 0.0543,-0.0511, 0.0075 }, {  5.039,-9.2324, 4.6003,-0.6253 },
    {-13.782, 36.397,-20.534, 2.9159 }, { 14.661,-42.962, 27.731,-4.1101 } 
  },
  // -------- Initial state nucleon-nucleon, outgoing h,K,Y --------
  { { 0.0716, 0.0926,-0.0515, 0.0058 }, {  3.096,-3.2186, 0.8989,-0.0017 },
    {-11.125, 20.273,-7.5084, 0.7022 }, {  18.13,-33.245, 13.188,-1.4854 } 
  },
  // -------- Initial state (h,K,Y,g)-nucleon, outgoing N --------
  { { 0.1729, -0.145, 0.0454,-0.0048 }, {  7.108,-13.032, 8.3515,-1.4095 },
    {-17.961, 41.781, -30.26, 5.3505 }, { 16.403,-40.799, 32.882,-6.0946 } 
  },
  // -------- Initial state (h,K,Y,g)-nucleon, outgoing h,K,Y --------
  { { 0.0376, 0.2383,-0.1541,  0.025 }, { 1.4331, 1.8253,-1.5201, 0.3059 },
    { -3.135, 1.7648,-1.5692, 0.3252 }, { 6.4864,-16.735, 17.185,-3.5277 } 
  }
};
