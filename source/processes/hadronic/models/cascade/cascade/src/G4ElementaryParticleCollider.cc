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

#include "G4ElementaryParticleCollider.hh"
#include "G4CascadeChannel.hh"
#include "G4CascadeChannelTables.hh"
#include "G4CascadeInterpolator.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclParticleNames.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4ParticleLargerEkin.hh"
#include "Randomize.hh"
#include <algorithm>
#include <cfloat>
#include <vector>

using namespace G4InuclParticleNames;
using namespace G4InuclSpecialFunctions;

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;


G4ElementaryParticleCollider::G4ElementaryParticleCollider()
  : G4CascadeColliderBase("G4ElementaryParticleCollider") {}


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
    printFinalStateTables();		// For diagnostic reporting
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

      G4LorentzVector mom = sampleCMmomentumFor2to2(is, kw, ekin, pmod);

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
G4ElementaryParticleCollider::getMomModuleFor2toMany(G4int is, G4int /*mult*/, 
					             G4int knd, 
					     	     G4double ekin) const 
{
  if (verboseLevel > 2) {
    G4cout << " >>> G4ElementaryParticleCollider::getMomModuleFor2toMany "
	   << " is " << is << " knd " << knd << " ekin " << ekin << G4endl;
  }

  G4double S = inuclRndm();
  G4double PS = 0.0;
  G4double PR = 0.0;
  G4double PQ = 0.0;
  G4int KM = 2;
  G4int IL = 4;
  G4int JK = 4;
  G4int JM = 2;
  G4int IM = 3;

  if (is == 1 || is == 2 || is == 4) KM = 1;
  if (knd == 1 || knd == 2) JK = 0;

  if (verboseLevel > 3) {
    G4cout << " S " << S << " KM " << KM << " IL " << IL << " JK " << JK
	   << " JM " << JM << " IM " << IM << G4endl;
  }

  for(G4int i = 0; i < 4; i++) {
    G4double V = 0.0;
    for(G4int k = 0; k < 4; k++) {
      if (verboseLevel > 3) {
	G4cout << " k " << k << " : rmn[k+JK][i+IL][KM-1] "
	       << rmn[k+JK][i+IL][KM-1] << " ekin^k " << std::pow(ekin, k)
	       << G4endl;
      }

      V += rmn[k + JK][i + IL][KM - 1] * std::pow(ekin, k);
    }

    if (verboseLevel > 3) {
      G4cout << " i " << i << " : V " << V << " S^i " << std::pow(S, i)
	     << G4endl;
    }

    PR += V * std::pow(S, i);
    PQ += V;
  }

  if (verboseLevel > 3) G4cout << " PR " << PR << " PQ " << PQ << G4endl;

  if (knd == 1 || knd == 2) JM = 1;
  if (verboseLevel > 3) G4cout << " JM " << JM << G4endl;

  for(G4int im = 0; im < 3; im++) {
    if (verboseLevel >3) {
      G4cout << " im " << im << " : rmn[8+IM+im][7+JM][KM-1] "
	     << rmn[8+IM+im][7+JM][KM-1] << " ekin^im " << std::pow(ekin, im)
	     << G4endl;
    }
    PS += rmn[8 + IM + im][7 + JM][KM - 1] * std::pow(ekin, im);
  }
  
  G4double PRA = PS * std::sqrt(S) * (PR + (1 - PQ) * (S*S*S*S));

  if (verboseLevel > 3) 
    G4cout << " PS " << PS << " PRA = PS*sqrt(S)*(PR+(1-PQ)*S^4) " << PRA
	   << G4endl;

  return std::fabs(PRA);
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
  G4int K = 3;
  G4int J = 1;

  if(is == 1 || is == 2 || is == 4) K = 1;

  if(knd == 1 || knd == 2) J = 0;

  G4int itry = 0;

  while(std::fabs(ct) > 1.0 && itry < itry_max) {
    itry++;
    G4double S = inuclRndm();
    G4double U = 0.0;
    G4double W = 0.0;

    for(G4int l = 0; l < 4; l++) {
      G4double V = 0.0;

      for(G4int im = 0; im < 4; im++) {
	V += abn[im][l][K+J-1] * std::pow(ekin, im);
      };

      U += V;
      W += V * std::pow(S, l);
    };  
    ct = 2.0 * std::sqrt(S) * (W + (1.0 - U) * (S*S*S*S)) - 1.0;
  };

  if(itry == itry_max) {

    if(verboseLevel > 2){
      G4cout << " particleSCMmomentumFor2to3 -> itry = itry_max " << itry << G4endl;
    }

    ct = 2.0 * inuclRndm() - 1.0;
  };

  return generateWithFixedTheta(ct, pmod);
}


G4LorentzVector 
G4ElementaryParticleCollider::sampleCMmomentumFor2to2(G4int is, G4int kw, 
                                                      G4double ekin, 
			                              G4double pscm) const 
{
  if (verboseLevel > 3)
    G4cout << " >>> G4ElementaryParticleCollider::sampleCMmomentumFor2to2" 
	   << " is " << is << " kw " << kw << " ekin " << ekin << " pscm "
	   << pscm << G4endl;

  G4double pA=0.0, pC=0.0, pCos=0.0, pFrac=0.0;		// Angular parameters

  // Arrays below are parameters for two-exponential sampling of angular
  // distributions of two-body scattering in the specified channels

  if (is == 1 || is == 2 || is == 4 ||
      is == 21 || is == 23 || is == 25 || is == 27 || is ==29 || is == 31 ||
      is == 42 || is == 46 || is == 50 || is == 54 || is ==58 || is == 62) {
    // nucleon-nucleon or hyperon-nucleon
    if (verboseLevel > 3) G4cout << " nucleon/hyperon elastic" << G4endl;

    static const G4double nnke[9] =  { 0.0,   0.44, 0.59,   0.80,   1.00,   2.24,   4.40,   6.15,  10.00};
    static const G4double nnA[9] =   { 0.0,   0.34, 2.51,   4.59,   4.2,    5.61,   6.38,   7.93,   8.7};
    static const G4double nnC[9] =   { 0.0,   0.0,  1.21,   1.54,   1.88,   1.24,   1.91,   4.04,   8.7};
    static const G4double nnCos[9] = {-1.0,  -1.0, 0.4226, 0.4226, 0.4384, 0.7193, 0.8788, 0.9164,  0.95};
    static const G4double nnFrac[9] = {1.0,   1.0, 0.4898, 0.7243, 0.7990, 0.8892, 0.8493, 0.9583,  1.0};

    static G4CascadeInterpolator<9> interp(nnke);	// Only need one!
    pA = interp.interpolate(ekin, nnA);
    pC = interp.interpolate(ekin, nnC);
    pCos = interp.interpolate(ekin, nnCos);
    pFrac = interp.interpolate(ekin, nnFrac);

  } else if (kw == 2 && (is == 9 || is == 18)) {
    // gamma p -> pi+ n, gamma p -> pi0 p, gamma p -> K Y (and isospin variants)
    // for now and due to lack of better data, use the gamma p -> pi+ n angular
    // distribution for all of these channels
    if (verboseLevel > 3)
      G4cout << " gamma-nucleon inelastic with 2-body final state" << G4endl;

    static const G4double gnke[10] =   {0.0,   0.11,  0.22,   0.26,  0.30,  0.34,  0.42,   0.59,   0.79,  10.0};
    static const G4double gnA[10] =    {0.0,   0.0,   5.16,   5.55,  5.33,  7.40, 13.55,  13.44,  13.31,   7.3};
    static const G4double gnC[10] =    {0.0, -10.33, -5.44,  -5.92, -4.27, -0.66,  1.37,   1.07,   0.52,   7.3};
    static const G4double gnCos[10] =  {1.0,   1.0,   0.906,  0.940, 0.940, 0.906, 0.906,  0.91,   0.91,   0.94};
    static const G4double gnFrac[10] = {0.0,   0.0,   0.028,  0.012, 0.014, 0.044, 0.087,  0.122,  0.16,   1.0};

    static G4CascadeInterpolator<10> interp(gnke);
    pA = interp.interpolate(ekin, gnA);
    pC = interp.interpolate(ekin, gnC);
    pCos = interp.interpolate(ekin, gnCos);
    pFrac = interp.interpolate(ekin, gnFrac);

  } else if (kw == 2) {
    // pi- p -> pi0 n, pi0 p -> pi+ n, pi- p -> K Y, pi0 p -> K Y (and isospin variants)
    // includes charge and strangeness exchange  
    if (verboseLevel > 3)
      G4cout << " pion-nucleon inelastic with 2-body final state " << G4endl;

    static const G4double qxke[10] =   {0.0,   0.062,  0.12,   0.217,  0.533,  0.873,  1.34,   2.86,   5.86,  10.0};
    static const G4double qxA[10] =    {0.0,   0.0,    2.48,   7.93,  10.0,    9.78,   5.08,   8.13,   8.13,   8.13};
    static const G4double qxC[10] =    {0.0, -39.58, -12.55,  -4.38,   1.81,  -1.99,  -0.33,   1.2,    1.43,   8.13};
    static const G4double qxCos[10] =  {1.0,   1.0,    0.604, -0.033,  0.25,   0.55,   0.65,   0.80,   0.916,  0.916};
    static const G4double qxFrac[10] = {0.0,   0.0,    0.1156, 0.5832, 0.8125, 0.3357, 0.3269, 0.7765, 0.8633, 1.0};

    static G4CascadeInterpolator<10> interp(qxke);	// Only need one!
    pA = interp.interpolate(ekin, qxA);
    pC = interp.interpolate(ekin, qxC);
    pCos = interp.interpolate(ekin, qxCos);
    pFrac = interp.interpolate(ekin, qxFrac);

  } else if (is == 3 || is == 7 || is == 9 || is == 11 || is == 17 ||
             is == 10 || is == 14 || is == 18 || is == 26 || is == 30) {
    // pi+p, pi0p, gammap, k+p, k0bp, pi-n, pi0n, gamman, k-n, or k0n
    if (verboseLevel > 3) G4cout << " meson-nucleon elastic (1)" << G4endl;

    static const G4double hn1ke[10] =   {0.0,  0.062,  0.12,   0.217,  0.533,  0.873,  1.34,   2.86,   5.86,  10.0};
    static const G4double hn1A[10] =    {0.0,  0.0,   27.58,  19.83,   6.46,   4.59,   6.47,   6.68,   6.43,   6.7};
    static const G4double hn1C[10] =    {0.0, -26.4, -30.55, -19.42,  -5.05,  -5.24,  -1.00,   2.14,   2.9,    6.7};
    static const G4double hn1Cos[10] =  {1.0,  1.0,    0.174, -0.174, -0.7,   -0.295,  0.5,    0.732,  0.837,  0.89};
    static const G4double hn1Frac[10] = {0.0,  0.0,    0.2980, 0.7196, 0.9812, 0.8363, 0.5602, 0.9601, 0.9901, 1.0};

    static G4CascadeInterpolator<10> interp(hn1ke);	// Only need one!
    pA = interp.interpolate(ekin, hn1A);
    pC = interp.interpolate(ekin, hn1C);
    pCos = interp.interpolate(ekin, hn1Cos);
    pFrac = interp.interpolate(ekin, hn1Frac);

  } else if (is == 5 || is == 6 || is == 13 || is == 34 || is == 22 ||
	     is == 15) {
    // pi-p, pi+n, k-p, k0bn, k+n, or k0p
    if (verboseLevel > 3) G4cout << " meson-nucleon elastic (2)" << G4endl;

    static const G4double hn2ke[10] =   {0.0,  0.062, 0.12,   0.217,  0.533,  0.873,  1.34,   2.86,   5.86,  10.0};
    static const G4double hn2A[10] =    {0.0, 27.08, 19.32,   9.92,   7.74,   9.86,   5.51,   7.25,   7.23,   7.3};
    static const G4double hn2C[10] =    {0.0,  0.0, -19.49, -15.78,  -9.78,  -2.74,  -1.16,   2.31,   2.96,   7.3};
    static const G4double hn2Cos[10] = {-1.0, -1.0,  -0.235, -0.259, -0.276,  0.336,  0.250,  0.732,  0.875,  0.9};
    static const G4double hn2Frac[10] = {1.0,  1.0,   0.6918, 0.6419, 0.7821, 0.6542, 0.8382, 0.9722, 0.9784, 1.0};

    static G4CascadeInterpolator<10> interp(hn2ke);	// Only need one!
    pA = interp.interpolate(ekin, hn2A);
    pC = interp.interpolate(ekin, hn2C);
    pCos = interp.interpolate(ekin, hn2Cos);
    pFrac = interp.interpolate(ekin, hn2Frac);

  } else {
    if (verboseLevel)
      G4cerr << " G4ElementaryParticleCollider::sampleCMmomentumFor2to2:"
	     << " interaction is=" << is << " not recognized " << G4endl;
  } 

  // Bound parameters by their physical ranges
  pCos = std::max(-1.,std::min(pCos,1.));
  pFrac = std::max(0.,std::min(pFrac,1.));

  // Use parameters determined above to get polar angle
  G4double ct = sampleCMcosFor2to2(pscm, pFrac, pA, pC, pCos);

  return generateWithFixedTheta(ct, pscm);
}


G4double
G4ElementaryParticleCollider::sampleCMcosFor2to2(G4double pscm, G4double pFrac,
                                                 G4double pA, G4double pC,
                                                 G4double pCos) const 
{
  if (verboseLevel>3) {
    G4cout << " sampleCMcosFor2to2: pscm " << pscm << " pFrac " << pFrac
	   << " pA " << pA << " pC " << pC << " pCos " << pCos << G4endl;
  }

  G4bool smallAngle = (G4UniformRand() < pFrac);	// 0 < theta < theta0

  G4double term1 = 2.0 * pscm*pscm * (smallAngle ? pA : pC);

  if (std::abs(term1) < 1e-7) return 1.0;	// No actual scattering here!
  if (term1 > DBL_MAX_EXP) return 1.0;

  G4double term2 = std::exp(-2.0*term1);
  G4double randScale = (std::exp(-term1*(1.0 - pCos)) - term2)/(1.0 - term2);

  G4double randVal;
  if (smallAngle) randVal = (1.0 - randScale)*G4UniformRand() + randScale;
  else randVal = randScale*G4UniformRand();

  G4double costheta = 1.0 + std::log(randVal*(1.0 - term2) + term2)/term1;

  if (verboseLevel>3) {
    G4cout << " term1 " << term1 << " term2 " << term2 << " randVal "
	   << randVal << " => costheta " << costheta << G4endl;
  }

  return costheta;
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


// Dump lookup tables for N-body final states

void G4ElementaryParticleCollider::
printFinalStateTables(std::ostream& os) const {
  G4CascadeChannelTables::PrintTable(pro*pro, os);
  G4CascadeChannelTables::PrintTable(neu*pro, os);
  G4CascadeChannelTables::PrintTable(neu*neu, os);
  G4CascadeChannelTables::PrintTable(kmi*neu, os);
  G4CascadeChannelTables::PrintTable(kmi*pro, os);
  G4CascadeChannelTables::PrintTable(kpl*neu, os);
  G4CascadeChannelTables::PrintTable(kpl*pro, os);
  G4CascadeChannelTables::PrintTable(k0b*neu, os);
  G4CascadeChannelTables::PrintTable(k0b*pro, os);
  G4CascadeChannelTables::PrintTable(k0*neu, os);
  G4CascadeChannelTables::PrintTable(k0*pro, os);
  G4CascadeChannelTables::PrintTable(lam*neu, os);
  G4CascadeChannelTables::PrintTable(lam*pro, os);
  G4CascadeChannelTables::PrintTable(pim*neu, os);
  G4CascadeChannelTables::PrintTable(pim*pro, os);
  G4CascadeChannelTables::PrintTable(pip*neu, os);
  G4CascadeChannelTables::PrintTable(pip*pro, os);
  G4CascadeChannelTables::PrintTable(pi0*neu, os);
  G4CascadeChannelTables::PrintTable(pi0*pro, os);
  G4CascadeChannelTables::PrintTable(sm*neu, os);
  G4CascadeChannelTables::PrintTable(sm*pro, os);
  G4CascadeChannelTables::PrintTable(sp*neu, os);
  G4CascadeChannelTables::PrintTable(sp*pro, os);
  G4CascadeChannelTables::PrintTable(s0*neu, os);
  G4CascadeChannelTables::PrintTable(s0*pro, os);
  G4CascadeChannelTables::PrintTable(xim*neu, os);
  G4CascadeChannelTables::PrintTable(xim*pro, os);
  G4CascadeChannelTables::PrintTable(xi0*neu, os);
  G4CascadeChannelTables::PrintTable(xi0*pro, os);

  os << " * * * PRELIMINARY -- GAMMA-NUCLEON TABLES * * *" << G4endl;
  G4CascadeChannelTables::PrintTable(gam*neu, os);
  G4CascadeChannelTables::PrintTable(gam*pro, os);
}


// Parameter array for momentum calculation of many body final states
const G4double G4ElementaryParticleCollider::rmn[14][10][2] = {
  {{0.5028,   0.6305}, {3.1442, -3.7333}, {-7.8172,  13.464}, {8.1667, -18.594}, 
   {1.6208,   1.9439}, {-4.3139, -4.6268}, {12.291,  9.7879}, {-15.288, -9.6074}, 
   {   0.0,     0.0}, {   0.0,      0.0}},     

  {{0.9348,   2.1801}, {-10.59,  1.5163}, { 29.227,  -16.38}, {-34.55,  27.944}, 
   {-0.2009, -0.3464}, {1.3641,   1.1093}, {-3.403, -1.9313}, { 3.8559,  1.7064}, 
   {   0.0,     0.0}, {    0.0,     0.0}},    
  
  {{-0.0967, -1.2886}, {4.7335,  -2.457}, {-14.298,  15.129}, {17.685, -23.295}, 
   { 0.0126,  0.0271}, {-0.0835, -0.1164}, { 0.186,  0.2697}, {-0.2004, -0.3185}, 
   {   0.0,     0.0}, {    0.0,     0.0}},    
  
  {{-0.025,   0.2091}, {-0.6248, 0.5228}, { 2.0282, -2.8687}, {-2.5895, 4.2688}, 
   {-0.0002, -0.0007}, {0.0014,   0.0051}, {-0.0024, -0.015}, { 0.0022,  0.0196}, 
   {    0.0,    0.0}, {    0.0,     0.0}},     
  
  {{1.1965,   0.9336}, {-0.8289,-1.8181}, { 1.0426,  5.5157}, { -1.909,-8.5216}, 
   { 1.2419,  1.8693}, {-4.3633, -5.5678}, { 13.743, 14.795}, {-18.592, -16.903}, 
   {    0.0,    0.0}, {    0.0,     0.0}},     
  
  {{0.287,    1.7811}, {-4.9065,-8.2927}, { 16.264,  20.607}, {-19.904,-20.827}, 
   {-0.244,  -0.4996}, {1.3158,   1.7874}, {-3.5691, -4.133}, { 4.3867,  3.8393}, 
   {    0.0,    0.0}, {   0.0,      0.0}}, 
  
  {{-0.2449, -1.5264}, { 2.9191, 6.8433}, {-9.5776, -16.067}, { 11.938, 16.845}, 
   {0.0157,   0.0462}, {-0.0826, -0.1854}, { 0.2143, 0.4531}, {-0.2585, -0.4627}, 
   {    0.0,    0.0}, {   0.0,      0.0}},
  
  {{0.0373,   0.2713}, {-0.422, -1.1944}, { 1.3883,  2.7495}, {-1.7476,-2.9045}, 
   {-0.0003, -0.0013}, {0.0014,   0.0058}, {-0.0034,-0.0146}, { 0.0039,  0.0156}, 
   {    0.0,    0.0}, {    0.0,     0.0}},     
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, 
   { 0.1451,  0.0929},{ 0.1538,  0.1303}},  
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, 
   { 0.4652,  0.5389},{ 0.2744,  0.4071}},  
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0},
   { -0.033, -0.0545},{-0.0146, -0.0288}},  
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0},
   { 0.6296,  0.1491},{ 0.8381,  0.1802}},  
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0},
   { 0.1787,   0.385},{ 0.0086,  0.3302}},  
  
  {{   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0}, {    0.0,     0.0},
   {    0.0,     0.0}, {   0.0,      0.0}, {    0.0,    0.0}, {    0.0,     0.0},
   {-0.0026, -0.0128},{ 0.0033, -0.0094}}
};

const G4double G4ElementaryParticleCollider::abn[4][4][4] = {
  {{0.0856,  0.0716,  0.1729,  0.0376},  {5.0390,  3.0960,  7.1080,  1.4331},
   {-13.782, -11.125, -17.961, -3.1350},  {14.661,  18.130,  16.403,  6.4864}},
  {{0.0543,  0.0926, -0.1450,  0.2383}, {-9.2324, -3.2186, -13.032,  1.8253},
   {36.397,  20.273,  41.781,  1.7648}, {-42.962, -33.245, -40.799, -16.735}},
  {{-0.0511, -0.0515,  0.0454, -0.1541}, {4.6003,  0.8989,  8.3515, -1.5201},
   {-20.534, -7.5084, -30.260, -1.5692},  {27.731,  13.188,  32.882,  17.185}},
  {{0.0075,  0.0058, -0.0048,  0.0250}, {-0.6253, -0.0017, -1.4095,  0.3059},
   {2.9159,  0.7022,  5.3505,  0.3252}, {-4.1101, -1.4854, -6.0946, -3.5277}} 
};
