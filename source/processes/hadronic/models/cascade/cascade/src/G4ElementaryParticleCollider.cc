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

#include "G4CascadeKplusPChannel.hh"
#include "G4CascadeKplusNChannel.hh"
#include "G4CascadeKzeroPChannel.hh"
#include "G4CascadeKzeroNChannel.hh"
#include "G4CascadeKminusPChannel.hh"
#include "G4CascadeKminusNChannel.hh"
#include "G4CascadeKzeroBarPChannel.hh"
#include "G4CascadeKzeroBarNChannel.hh"
#include "G4CascadeLambdaPChannel.hh"
#include "G4CascadeLambdaNChannel.hh"
#include "G4CascadeSigmaPlusPChannel.hh"
#include "G4CascadeSigmaPlusNChannel.hh"
#include "G4CascadeSigmaZeroPChannel.hh"
#include "G4CascadeSigmaZeroNChannel.hh"
#include "G4CascadeSigmaMinusPChannel.hh"
#include "G4CascadeSigmaMinusNChannel.hh"
#include "G4CascadeXiZeroPChannel.hh"
#include "G4CascadeXiZeroNChannel.hh"
#include "G4CascadeXiMinusPChannel.hh"
#include "G4CascadeXiMinusNChannel.hh"

#include "G4Collider.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4ParticleLargerEkin.hh"
#include <algorithm>

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;


G4ElementaryParticleCollider::G4ElementaryParticleCollider()
  : verboseLevel(1)
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider ctor " << G4endl;
  }

  initializeArrays();
}


G4CollisionOutput  
G4ElementaryParticleCollider::collide(G4InuclParticle* bullet,
				      G4InuclParticle* target) 
{
  G4InuclElementaryParticle* particle1 =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* particle2 =	
    dynamic_cast<G4InuclElementaryParticle*>(target);

  G4CollisionOutput output;

  if (!(particle1 && particle2)) {
    G4cout << " ElementaryParticleCollider -> can collide only particle with particle " 
           << G4endl;
  } else {
    collide(particle1, particle2, output);
  }
  return output;	
}


void 
G4ElementaryParticleCollider::collide(G4InuclElementaryParticle* particle1,
				      G4InuclElementaryParticle* particle2,
				      G4CollisionOutput& output) {

  // Generate nucleon or pion collision with nucleon
  // or pion with quasi-deuteron

  if (!particle1->photon() && !particle2->photon()) { // ok
    if (particle1->nucleon() || particle2->nucleon()) { // ok
      G4LorentzConvertor convertToSCM;
      if(particle2->nucleon()) {
	convertToSCM.setBullet(particle1->getMomentum(), particle1->getMass());
	convertToSCM.setTarget(particle2->getMomentum(), particle2->getMass());
      } else {
	convertToSCM.setBullet(particle2->getMomentum(), particle2->getMass());
	convertToSCM.setTarget(particle1->getMomentum(), particle1->getMass());
      };  
      convertToSCM.toTheCenterOfMass();
      G4double ekin = convertToSCM.getKinEnergyInTheTRS();
      G4double etot_scm = convertToSCM.getTotalSCMEnergy();
      G4double pscm = convertToSCM.getSCMMomentum();

      std::vector<G4InuclElementaryParticle> particles = 	    
	generateSCMfinalState(ekin, etot_scm, pscm, particle1, particle2, &convertToSCM);

      if(verboseLevel > 2){
	G4cout << " particles " << particles.size() << G4endl;

	for(G4int i = 0; i < G4int(particles.size()); i++) 
	  particles[i].printParticle();

      }
      if(!particles.empty()) { // convert back to Lab
	particleIterator ipart;
	for(ipart = particles.begin(); ipart != particles.end(); ipart++) {	
	  G4CascadeMomentum mom = 
	    convertToSCM.backToTheLab(ipart->getMomentum());
	  ipart->setMomentum(mom); 
	};
	std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
      };
      output.addOutgoingParticles(particles);

    } else {
      if(particle1->quasi_deutron() || particle2->quasi_deutron()) {
	if(particle1->pion() || particle2->pion()) {
	  G4LorentzConvertor convertToSCM;
	  if(particle1->pion()) {
	    convertToSCM.setBullet(particle1->getMomentum(), particle1->getMass());
	    convertToSCM.setTarget(particle2->getMomentum(), particle2->getMass());
	  } else {
	    convertToSCM.setBullet(particle2->getMomentum(), particle2->getMass());
	    convertToSCM.setTarget(particle1->getMomentum(), particle1->getMass());
	  }; 
	  convertToSCM.toTheCenterOfMass(); 
	  G4double etot_scm = convertToSCM.getTotalSCMEnergy();
	  std::vector<G4InuclElementaryParticle> particles = 
	    generateSCMpionAbsorption(etot_scm, particle1, particle2);

	  if(!particles.empty()) { // convert back to Lab
	    particleIterator ipart;
	    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	      G4CascadeMomentum mom = 
		convertToSCM.backToTheLab(ipart->getMomentum());
	      ipart->setMomentum(mom); 
	    };
	    std::sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
	    output.addOutgoingParticles(particles);
	  };

	} else {
	  G4cout << " ElementaryParticleCollider -> can only collide pions with deuterons " 
                 << G4endl;
	};
      } else {
	G4cout << " ElementaryParticleCollider -> can only collide something with nucleon or deuteron " 
               << G4endl;
      };
    };  

  } else {

    G4cout << " ElementaryParticleCollider -> cannot collide photons " 
           << G4endl;
  }; 

}


G4int 
G4ElementaryParticleCollider::generateMultiplicity(G4int is, 
						   G4double ekin) const 
{
  G4int mul = 0;
  G4int l = is;

  if ( ( l > 10 && l < 14 ) || ( l > 14 && l < 63 ) ) {
    // strange particle branch
    if ( l == 11 ) {
      mul = G4CascadeKplusPChannel::getMultiplicity(ekin);
    } else if ( l == 13 ) {
      mul = G4CascadeKminusPChannel::getMultiplicity(ekin);
    } else if ( l == 15 ) {
      mul = G4CascadeKzeroPChannel::getMultiplicity(ekin);
    } else if ( l == 17 ) {
      mul = G4CascadeKzeroBarPChannel::getMultiplicity(ekin);
    } else if ( l == 21 ) {
      mul = G4CascadeLambdaPChannel::getMultiplicity(ekin);
    } else if ( l == 23 ) {
      mul = G4CascadeSigmaPlusPChannel::getMultiplicity(ekin);
    } else if ( l == 25 ) {
      mul = G4CascadeSigmaZeroPChannel::getMultiplicity(ekin);
    } else if ( l == 27 ) {
      mul = G4CascadeSigmaMinusPChannel::getMultiplicity(ekin);
    } else if ( l == 29 ) {
      mul = G4CascadeXiZeroPChannel::getMultiplicity(ekin);
    } else if ( l == 31 ) {
      mul = G4CascadeXiMinusPChannel::getMultiplicity(ekin);

    } else if ( l == 22 ) {
      mul = G4CascadeKplusNChannel::getMultiplicity(ekin);
    } else if ( l == 26 ) {
      mul = G4CascadeKminusNChannel::getMultiplicity(ekin);
    } else if ( l == 30 ) {
      mul = G4CascadeKzeroNChannel::getMultiplicity(ekin);
    } else if ( l == 34 ) {
      mul = G4CascadeKzeroBarNChannel::getMultiplicity(ekin);
    } else if ( l == 42 ) {
      mul = G4CascadeLambdaNChannel::getMultiplicity(ekin);
    } else if ( l == 46 ) {
      mul = G4CascadeSigmaPlusNChannel::getMultiplicity(ekin);
    } else if ( l == 50 ) {
      mul = G4CascadeSigmaZeroNChannel::getMultiplicity(ekin);
    } else if ( l == 54 ) {
      mul = G4CascadeSigmaMinusNChannel::getMultiplicity(ekin);
    } else if ( l == 58 ) {
      mul = G4CascadeXiZeroNChannel::getMultiplicity(ekin);
    } else if ( l == 62 ) {
      mul = G4CascadeXiMinusNChannel::getMultiplicity(ekin);

    } else {
      G4cout << " G4ElementaryParticleCollider:" 
             << " Unknown strange interaction channel - multiplicity not generated " 
             << G4endl;
    }

  // Non-strange branch 
  } else if (l == 1 || l == 4) {
    mul = nucSampler.GetMultiplicityT1(ekin) - 2;
    if (mul > 7) G4cout << " Nuc sampler pp mult too high: mul = " << mul << G4endl;
  } else if (l == 2) {
    mul = nucSampler.GetMultiplicityT0(ekin) - 2;
    if (mul > 7) G4cout << " Nuc sampler np mult too high: mul = " << mul << G4endl;
  } else if (l == 3 || l == 10) {
    // |T,Tz> = |3/2,3/2> 
    mul = piSampler.GetMultiplicityT33(ekin) - 2;
  } else if (l == 5 || l == 6) {
    // |T,Tz> = |3/2,1/2> 
    mul = piSampler.GetMultiplicityT31(ekin) - 2;
  } else if (l == 7 || l == 14) {
    // |T,Tz> = |1/2,1/2> 
    mul = piSampler.GetMultiplicityT11(ekin) - 2;
  } else {
    G4cout << " G4ElementaryParticleCollider: "
           << " Unknown interaction channel - multiplicity not generated " 
           << G4endl;
  }

  if(verboseLevel > 3){
    G4cout << " G4ElementaryParticleCollider::generateMultiplicity: "  
           << " multiplicity = " << mul + 2 << G4endl; 
  }

  return mul + 2;
}

 
std::vector<G4InuclElementaryParticle> 
G4ElementaryParticleCollider::generateSCMfinalState(G4double ekin, 
		                     G4double etot_scm, 
		                     G4double pscm,
		                     G4InuclElementaryParticle* particle1,
		                     G4InuclElementaryParticle* particle2, 
	                             G4LorentzConvertor* toSCM) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMfinalState" 
           << G4endl;
  }

  const G4double ang_cut = 0.9999;
  const G4double difr_const = 0.3678794;   
  const G4int itry_max = 10;
  G4InuclElementaryParticle dummy;
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4int> particle_kinds;
  G4int type1 = particle1->type();
  G4int type2 = particle2->type();
  G4int is = type1 * type2;

  if(verboseLevel > 3){
    G4cout << " is " << is << G4endl;
  }

  G4int multiplicity = 0;
  G4bool generate = true;
   
  while (generate) {

    if(multiplicity == 0) {
      multiplicity = generateMultiplicity(is, ekin);
    } else {
      multiplicity = generateMultiplicity(is, ekin);
      particle_kinds.resize(0);
    }

    if(multiplicity == 2) { // 2 -> 2
      G4int kw = 1;
      if ( (is > 10 && is < 14) || (is > 14 && is < 63) ) {
        particle_kinds =
	  generateStrangeChannelPartTypes(is, 2, ekin);

        G4int finaltype = particle_kinds[0]*particle_kinds[1];
        if (finaltype != is) kw = 2;  // Charge or strangeness exchange

      } else if (is == 1 || is == 2 || is == 4) {
          particle_kinds.push_back(type1);
          particle_kinds.push_back(type2);

      } else if (is == 3) {
        particle_kinds = piSampler.GetFSPartTypesForPipP(2, ekin);
        if (particle_kinds[0] != G4PionSampler::pro) kw = 2;

      } else if (is == 10) {
        particle_kinds = piSampler.GetFSPartTypesForPimN(2, ekin);
        if (particle_kinds[0] != G4PionSampler::neu) kw = 2;

      } else if (is == 5) {
        particle_kinds = piSampler.GetFSPartTypesForPimP(2, ekin);
        if (particle_kinds[0] != G4PionSampler::pro) kw = 2;

      } else if (is == 6) {
        particle_kinds = piSampler.GetFSPartTypesForPipN(2, ekin);
        if (particle_kinds[0] != G4PionSampler::neu) kw = 2;

      } else if (is == 7) {
        particle_kinds = piSampler.GetFSPartTypesForPizP(2, ekin);
        if (particle_kinds[0] != G4PionSampler::pro) kw = 2;

      } else if (is == 14) {
        particle_kinds = piSampler.GetFSPartTypesForPizN(2, ekin);
        if (particle_kinds[0] != G4PionSampler::neu) kw = 2;

      } else {
        G4cout << " Unexpected interaction type (2->2) is = " << is << G4endl;
      }

      G4CascadeMomentum mom;

      if (kw == 2) { // need to rescale momentum
	G4double m1 = dummy.getParticleMass(particle_kinds[0]);
	m1 *= m1;
	G4double m2 = dummy.getParticleMass(particle_kinds[1]);
	m2 *= m2;	 
	G4double a = 0.5 * (etot_scm * etot_scm - m1 - m2);
	G4double em = a * a - m1 * m2;

        if (em > 0) { //  see if it is possible to rescale?
	  G4double np = std::sqrt( em / (m1 + m2 + 2.0 * a));
	  mom = particleSCMmomentumFor2to2(is, kw, ekin, np);
	} else { // rescaling not possible
	  mom = particleSCMmomentumFor2to2(is, kw, ekin, pscm); 
	}

      } else {
	mom = particleSCMmomentumFor2to2(is, kw, ekin, pscm);
      };

      //      G4cout << " Particle kinds = " << particle_kinds[0] << " , " << particle_kinds[1] << G4endl;

      if (verboseLevel > 3){
	G4cout << " before rotation px " << mom[1] << " py " << mom[2] <<
	  " pz " << mom[3] << G4endl;
      }

      mom = toSCM->rotate(mom); 

      if (verboseLevel > 3){
	G4cout << " after rotation px " << mom[1] << " py " << mom[2] <<
	  " pz " << mom[3] << G4endl;
      }
      G4CascadeMomentum mom1;

      for (G4int i = 1; i < 4; i++) mom1[i] = -mom[i];

      particles.push_back(G4InuclElementaryParticle(mom, particle_kinds[0], 3));
      // register modelId
      particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[1],3));
      generate = false;

    } else { // 2 -> many

      if ( (is > 10 && is < 14) || (is > 14 && is < 63) ) {
        particle_kinds =
          generateStrangeChannelPartTypes(is, multiplicity, ekin);

      } else if (is == 1) {
        particle_kinds = nucSampler.GetFSPartTypesForPP(multiplicity, ekin);

      } else if (is == 2) {
        particle_kinds = nucSampler.GetFSPartTypesForNP(multiplicity, ekin);

      } else if (is == 4) {
        particle_kinds = nucSampler.GetFSPartTypesForNN(multiplicity, ekin);

      } else if (is == 3) {
        particle_kinds = piSampler.GetFSPartTypesForPipP(multiplicity, ekin);

      } else if (is == 10) {
        particle_kinds = piSampler.GetFSPartTypesForPimN(multiplicity, ekin);

      } else if (is == 5) {
        particle_kinds = piSampler.GetFSPartTypesForPimP(multiplicity, ekin);

      } else if (is == 6) {
        particle_kinds = piSampler.GetFSPartTypesForPipN(multiplicity, ekin);

      } else if (is == 7) {
        particle_kinds = piSampler.GetFSPartTypesForPizP(multiplicity, ekin);

      } else if (is == 14) {
        particle_kinds = piSampler.GetFSPartTypesForPizN(multiplicity, ekin);

      } else {
        G4cout << " Unexpected interaction type is = " << is << G4endl;
      }

      //      G4cout << " Particle kinds = " ;
      //      for (G4int i = 0; i < multiplicity; i++) G4cout << particle_kinds[i] << " , " ;
      //       G4cout << G4endl;

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
	  G4cout << " itry in while " << itry << G4endl;
	}

	std::vector<G4double> modules = 
	  generateMomModules(particle_kinds, multiplicity, is, ekin, etot_scm);

	if (G4int(modules.size()) == multiplicity) {

	  if (multiplicity == 3) { 
	    G4CascadeMomentum mom3 = 
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

	      G4CascadeMomentum mom1 = generateWithFixedTheta(ct, modules[0]);

	      mom1 = toSCM->rotate(mom3, mom1);

	      G4CascadeMomentum mom2;

	      for(G4int i = 1; i < 4; i++) mom2[i] = - (mom3[i] + mom1[i]);

	      bad = false;
	      generate = false;

	      particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[0]));
	      particles.push_back(G4InuclElementaryParticle(mom2, particle_kinds[1]));
	      particles.push_back(G4InuclElementaryParticle(mom3, particle_kinds[2]));
	    };

	  } else { // multiplicity > 3

	    // generate first mult - 2 momentums
	    std::vector<G4CascadeMomentum> scm_momentums;
	    G4CascadeMomentum tot_mom;

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
		  G4cout << " s1 * alf * std::exp(-s1 / p0) " << s1 * alf * std::exp(-s1 / p0) 
			 << " s2 " << s2 << G4endl;
		}

		if(s1 * alf * std::exp(-s1 / p0) > s2) st = s1 / modules[i];

	      }; 

	      if(verboseLevel > 3){
		G4cout << " itry1 " << itry1 << " i " << i << " st " << st << G4endl;
	      }

	      if(itry1 == itry_max) {

		if(verboseLevel > 2){
		  G4cout << " high energy angles generation: itry1 " << itry1 << G4endl;
		}

		st = 0.5 * inuclRndm();
	      };

	      G4double ct = std::sqrt(1.0 - st * st);

	      if(inuclRndm() > 0.5) ct = -ct;

	      G4double pt = modules[i]*st;
	      G4double phi = randomPHI();

	      G4CascadeMomentum mom;

	      mom[1] = pt * std::cos(phi);
	      mom[2] = pt * std::sin(phi);
	      mom[3] = modules[i] * ct;

	      for(G4int i = 1; i < 4; i++) tot_mom[i] += mom[i];		 

	      scm_momentums.push_back(mom);
	    }; 

	    // handle last two
	    G4double tot_mod = std::sqrt(tot_mom[1] * tot_mom[1] + 
					 tot_mom[2] * tot_mom[2] + tot_mom[3] * tot_mom[3]); 
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

	      G4CascadeMomentum mom = 
		generateWithFixedTheta(ct, modules[multiplicity - 2]);

	      mom = toSCM->rotate(tot_mom, mom);
	      scm_momentums.push_back(mom);

	      // and the last one
	      G4CascadeMomentum mom1;

	      for (i = 1; i < 4; i++) mom1[i] = -mom[i] - tot_mom[i];

	      scm_momentums.push_back(mom1);  
	      bad = false;
	      generate = false;

	      if (verboseLevel > 2){
		G4cout << " ok for mult " << multiplicity << G4endl;
	      }

	      for (i = 0; i < multiplicity; i++) {
		particles.push_back(G4InuclElementaryParticle(
							      scm_momentums[i], particle_kinds[i]));
	      };
	    };
	  }; 
	};
      };

      if (itry == itry_max) {

	if (verboseLevel > 2){
	  G4cout << " cannot generate the distr. for mult " << multiplicity  <<
	    G4endl << " and set it to " << multiplicity - 1 << G4endl;
	}

      };
    };
  }; 

  if (verboseLevel > 3) {
    G4cout << " <<< G4ElementaryParticleCollider::generateSCMfinalState" << G4endl;
  }

  return particles;
}

std::vector<G4double> 
G4ElementaryParticleCollider::generateMomModules(
		   const std::vector<G4int>& kinds, 
		   G4int mult, 
		   G4int is, 
		   G4double ekin, 
		   G4double etot_cm) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateMomModules" 
           << G4endl;
  }

  if (verboseLevel > 2){
    G4cout << " mult " << mult << " is " << is << " ekin " << ekin << " etot_cm " <<
      etot_cm << G4endl;
  }

  const G4int itry_max = 10;
  const G4double small = 1.e-10;
  G4InuclElementaryParticle dummy;
  G4int itry = 0;

  std::vector<G4double> modules(mult);
  std::vector<G4double> masses2(mult);

  for (G4int i = 0; i < mult; i++) {
    G4double mass = dummy.getParticleMass(kinds[i]);
    masses2[i] = mass * mass;
  };

  G4double mass_last = std::sqrt(masses2[mult - 1]);

  if (verboseLevel > 3){
    G4cout << " knd_last " << kinds[mult - 1] << " mlast " 
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
	getMomModuleFor2toMany(is, mult, kinds[i], ekin);

      if (pmod < small) break;
      eleft -= std::sqrt(pmod * pmod + masses2[i]);

      if (verboseLevel > 3){
	G4cout << " kp " << kinds[i] << " pmod " << pmod << " mass2 " 
               << masses2[i] << G4endl;
	G4cout << " x1 " << eleft - mass_last << G4endl;
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
	  if(satisfyTriangle(modules)) {
	    return modules;
	  }

	} else {
	  return modules;
	}
      }
    }
  }

  modules.resize(0);
  return modules;    
}


G4bool 
G4ElementaryParticleCollider::satisfyTriangle(
			const std::vector<G4double>& modules) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::satisfyTriangle" 
           << G4endl;
  }

  G4bool good = true;
  if(modules.size() == 3) {

    if(std::fabs(modules[1] - modules[2]) > modules[0] || 
       modules[0] > modules[1] + modules[2] ||
       std::fabs(modules[0] - modules[2]) > modules[1] ||
       modules[1] > modules[0] + modules[2] ||
       std::fabs(modules[0] - modules[1]) > modules[2] ||
       modules[2] > modules[1] + modules[0]) good = false;

  }

  return good;
}


std::vector<G4int> 
G4ElementaryParticleCollider::generateStrangeChannelPartTypes(
                              G4int is, G4int mult, G4double ekin) const
{
  std::vector<G4int> kinds;

  if (is == 11) {
    kinds = G4CascadeKplusPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 13) {
    kinds = G4CascadeKminusPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 15) {
    kinds = G4CascadeKzeroPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 17) {
    kinds = G4CascadeKzeroBarPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 21) {
    kinds = G4CascadeLambdaPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 23) {
    kinds = G4CascadeSigmaPlusPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 25) {
    kinds = G4CascadeSigmaZeroPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 27) {
    kinds = G4CascadeSigmaMinusPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 29) {
    kinds = G4CascadeXiZeroPChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 31) {
    kinds = G4CascadeXiMinusPChannel::getOutgoingParticleTypes(mult, ekin);

  } else if (is == 22) {
    kinds = G4CascadeKplusNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 26) {
    kinds = G4CascadeKminusNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 30) {
    kinds = G4CascadeKzeroNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 34) {
    kinds = G4CascadeKzeroBarNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 42) {
    kinds = G4CascadeLambdaNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 46) {
    kinds = G4CascadeSigmaPlusNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 50) {
    kinds = G4CascadeSigmaZeroNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 54) {
    kinds = G4CascadeSigmaMinusNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 58) {
    kinds = G4CascadeXiZeroNChannel::getOutgoingParticleTypes(mult, ekin);
  } else if (is == 62) {
    kinds = G4CascadeXiMinusNChannel::getOutgoingParticleTypes(mult, ekin);

  } else {
    G4cout << " G4ElementaryParticleCollider:"
	   << " Unknown strange interaction channel - outgoing kinds not generated " 
           << G4endl;
  }

  return kinds;
}


G4double 
G4ElementaryParticleCollider::getMomModuleFor2toMany(G4int is, G4int mult, 
					             G4int knd, 
					     	     G4double ekin) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::getMomModuleFor2toMany" 
           << G4endl;
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

  if(is == 1 || is == 2 || is == 4) KM = 1;
  if(mult == 3) { IM = 0; IL = 0; };
  if(knd == 1 || knd == 2) JK = 0;

  for(G4int i = 0; i < 4; i++) {
    G4double V = 0.0;
    for(G4int k = 0; k < 4; k++) V += rmn[k + JK][i + IL][KM - 1] * std::pow(ekin, k);
    PR += V * std::pow(S, i);
    PQ += V;
  }

  if(knd == 1 || knd == 2) JM = 1;
  for(G4int m = 0; m < 3; m++) PS += rmn[8 + IM + m][7 + JM][KM - 1] * std::pow(ekin, m);
  G4double PRA = PS * std::sqrt(S) * (PR + (1 - PQ) * std::pow(S, 4));

  return std::fabs(PRA);
}


G4CascadeMomentum 
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

      for(G4int m = 0; m < 4; m++) {
	V += abn[m][l][K + J - 1] * std::pow(ekin, m);
      };

      U += V;
      W += V * std::pow(S, l);
    };  
    ct = 2.0 * std::sqrt(S) * (W + (1.0 - U) * std::pow(S, 4)) - 1.0;
  };

  if(itry == itry_max) {

    if(verboseLevel > 2){
      G4cout << " particleSCMmomentumFor2to3 -> itry = itry_max " << itry << G4endl;
    }

    ct = 2.0 * inuclRndm() - 1.0;
  };

  G4double pt = pmod * std::sqrt(1.0 - ct * ct);
  G4double phi = randomPHI();

  G4CascadeMomentum mom;

  mom[1] = pt * std::cos(phi);
  mom[2] = pt * std::sin(phi);
  mom[3] = pmod * ct;
  
  return mom;  
}


std::pair<G4double, G4double> 
G4ElementaryParticleCollider::adjustIntervalForElastic(G4double ekin, G4double ak, 
			                               G4double ae, G4int k, 
			                               G4int l, 
			                               const std::vector<G4double>& ssv,
			                               G4double st) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::adjustIntervalForElastic" 
           << G4endl;
  }

  const G4int itry_max = 100;
  const G4double small = 1.0e-4;
  G4double a = 1.0;
  G4double b = 0.0;
  G4int adj_type = 0;
  G4double s1 = 0.0;
  G4double s2 = 0.0;

  if (k == 1) {
    adj_type = 1;
    s1 = 0.0;
    s2 = 0.5;
  } else if(k == 2) {

    if(l != 2) {
      // adj_type == 0;

    } else { 

      adj_type = 1;
      s1 = 0.0;
      s2 = 0.67;
    };
        
  } else {

    if(ekin < 0.32) {
      // adj_type == 0;

    } else {

      adj_type = 2;
      s1 = 0.1813;
      s2 = 1.0;
    };     

  };

  if(adj_type > 0) {
    G4int itry = 0;
    G4double su;
    G4double ct;
  
    if(adj_type == 1) {
      G4double s2_old = s2;
      G4double s1c = s1;
      G4double s2_new;

      while(itry < itry_max) {
	itry++;
	s2_new = 0.5 * (s2_old + s1c);
	su = 0.0;      

	for(G4int i = 0; i < 4; i++) su += ssv[i] * std::pow(s2_new, i);

	ct = ak * std::sqrt(s2_new) * (su + (1.0 - st) * std::pow(s2_new, 4)) + ae;

	if(ct > 1.0) {
	  s2_old = s2_new;

	} else {

	  if(1.0 - ct < small) {

	    break;

	  } else {

	    s1c = s2_new;
	  };   
	}; 
      };

      a = s2_new - s1;
      b = s1;

    } else {

      G4double s1_old = s1;
      G4double s2c = s2;
      G4double s1_new = 0.0;

      while (itry < itry_max) {
	itry++;
	s1_new = 0.5 * (s1_old + s2c);
	su = 0.0;
      
	for (G4int i = 0; i < 4; i++) su += ssv[i] * std::pow(s1_new, i);
	ct = ak * std::sqrt(s1_new) * (su + (1.0 - st) * std::pow(s1_new, 4)) + ae;

	if(ct < -1.0) {
	  s1_old = s1_new;

	} else {

	  if (1.0 + ct < small) {

	    break;

	  } else {

	    s2c = s1_new;
	  };   
	}; 
      };
      a = s2 - s1_new;
      b = s1_new;
    }; 

    if (itry == itry_max) {

      if (verboseLevel > 2){
	G4cout << " in adjustIntervalForElastic: " << itry_max << G4endl;
	G4cout << " e " << ekin << " ak " << ak << " ae " << ae << G4endl << " k " << k
	       << " is " << l << " adj_type " << adj_type << G4endl; 
      }

      a = 1.0;
      b = 0.0;
    };
  };

  return std::pair<G4double, G4double>(a, b);
}  


G4CascadeMomentum 
G4ElementaryParticleCollider::particleSCMmomentumFor2to2(
			   G4int is, 
			   G4int kw, 
			   G4double ekin, 
			   G4double pscm) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::particleSCMmomentumFor2to2" 
           << G4endl;
  }

  const G4int itry_max = 100;
  const G4double ct_cut = 0.9999;
  const G4double huge_num = 60.0;
  G4int k = getElasticCase(is, kw, ekin);
  G4double ae = -1.0;
  G4double ak = 2.0;

  if(k == 1) {
    if(is != 2) { ak = 1.0; ae = 0.0;};
  } else if(k == 2) {
    if(is != 2) { ak = 0.5; ae = 0.0; };
  }

  G4double ct = 2.0;
  G4double ab;
  G4double ac;
  G4double ad;
  G4int itry = 0;
  
  if(k == 14) {
    ab = 7.5;
    if(is == 1 || is == 2 || is == 4) ab = 8.7;
    ac = -2.0 * ab * pscm * pscm;
    ad = 2.0 * ac;
    if(ad < -huge_num) {
      ad = std::exp(ad);

    } else {

      ad = std::exp(-huge_num);
    };   

    while(std::fabs(ct) > ct_cut && itry < itry_max) {
      itry++;
      ct = 1.0 - std::log(inuclRndm() * (1.0 - ad) + ad) / ac;
    };      

  } else if(k == 0) {
    ct = 2.0 * inuclRndm() - 1;

  } else { 
    G4int k1 = k - 1;
    // first set all coefficients
    std::vector<G4double> ssv(4);
    G4double st = 0.0;

    for(G4int i = 0; i < 4; i++) {
      G4double ss = 0.0;

      for(G4int m = 0; m < 4; m++) ss += ang[m][i][k1] * std::pow(ekin, m);
      st += ss;
      ssv[i] = ss;
    };

    G4double a = 1.0;
    G4double b = 0.0;

    if(k <= 3) {
      std::pair<G4double, G4double> ab = adjustIntervalForElastic(ekin, ak, ae, k, is, ssv, st);

      a = ab.first;
      b = ab.second;
    };

    while(std::fabs(ct) > ct_cut && itry < itry_max) {
      itry++;
      G4double mrand = a * inuclRndm() + b;

      G4double su = 0.0;

      for(G4int i = 0; i < 4; i++) su += ssv[i] * std::pow(mrand, i);

      ct = ak * std::sqrt(mrand) * (su + (1.0 - st) * std::pow(mrand, 4)) + ae;
    };

  };

  if(itry == itry_max) {
    if(verboseLevel > 2){
      G4cout << " particleSCMmomentumFor2to2 -> itry = itry_max " 
             << itry << G4endl;
    }
    ct = 2.0 * inuclRndm() - 1.0;
  }

  G4double pt = pscm * std::sqrt(1.0 - ct * ct);
  G4double phi = randomPHI();
  G4CascadeMomentum mom;

  mom[1] = pt * std::cos(phi);
  mom[2] = pt * std::sin(phi);
  mom[3] = pscm * ct;
  
  return mom;  
}


G4int 
G4ElementaryParticleCollider::getElasticCase(G4int is, 
					     G4int kw, 
					     G4double ekin) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::getElasticCase" 
           << G4endl;
  }

  G4int l = is;
  G4int k = 0; // isotropic

  if(l == 4) {
    l = 1;
  } else if(l == 10 || l == 7 || l == 14) {
    l = 3;
  } else if(l == 5 || l == 6) {
    l = 4;
  }

  if(l < 3) { // nucleon nucleon

    if(ekin > 2.8) {
      // DHW      k = 2;
      // DHW      if(ekin > 10.0) k = 14;
      k = 14;

    } else {
      if(l == 1) { // PP or NN
        if(ekin > 0.46) k = 1;
      } else {
        k = 3;
	if(ekin >= 0.97) k = 1;
      }
    }
  	 
  } else { // pi nucleon

    if(l == 3) { // pi+ P, pi- N, pi0 P, pi0 N
      k = 8;
      if(ekin > 0.08) k = 9;
      if(ekin > 0.3) k = 10;
      if(ekin > 1.0) k = 11;
      if(ekin > 2.4) k = 14;

    } else { // pi- P, pi+ N

      if(kw == 1) {
        k = 4;
        if(ekin > 0.08) k = 5;
        if(ekin > 0.3) k = 6;
        if(ekin > 1.0) k = 7;
        if (ekin > 2.4) k = 14;

      } else {

        k = 12;
	if (ekin > 0.08) k = 13;
        if (ekin > 0.3) k = 6;
        if (ekin > 1.0) k = 7;
        if (ekin > 2.4) k = 14;
      }
    }  
  }     

  return k;
}


std::vector<G4InuclElementaryParticle> 
G4ElementaryParticleCollider::generateSCMpionAbsorption(G4double etot_scm,
			             G4InuclElementaryParticle* particle1,
			             G4InuclElementaryParticle* particle2) const 
{
  if (verboseLevel > 3) {
    G4cout << " >>> G4ElementaryParticleCollider::generateSCMpionAbsorption" 
           << G4endl;
  }

  // generate nucleons momenta for pion absorption
  // the nucleon distribution assumed to be isotropic in SCM

  G4InuclElementaryParticle dummy;
  std::vector<G4InuclElementaryParticle> particles;
  std::vector<G4int> particle_kinds;
  G4int type1 = particle1->type();
  G4int type2 = particle2->type();

  // generate kinds

  if (type1 == 3) {

    if (type2 == 111) { // pi+ + PP -> ? 

      G4cout << " pion absorption: pi+ + PP -> ? " << G4endl;

      return particles;
    }
    else if(type2 == 112) { // pi+ + PN -> PP
      particle_kinds.push_back(1);
      particle_kinds.push_back(1);
    }
    else if(type2 == 122) { // pi+ + NN -> PN
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    };     
  }
  else if(type1 == 5) { 
    if(type2 == 111) { // pi- + PP -> PN 
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    }
    else if(type2 == 112) { // pi- + PN -> NN
      particle_kinds.push_back(2);
      particle_kinds.push_back(2);
    }
    else if(type2 == 122) { // pi- + NN -> ?

      G4cout << " pion absorption: pi- + NN -> ? " << G4endl;

      return particles;
    };     
  }
  else if(type1 == 7) {
    if(type2 == 111) { // pi0 + PP -> PP 
      particle_kinds.push_back(1);
      particle_kinds.push_back(1);
    }
    else if(type2 == 112) { // pi0 + PN -> PN
      particle_kinds.push_back(1);
      particle_kinds.push_back(2);
    }
    else if(type2 == 122) { // pi0 + NN -> ?
      particle_kinds.push_back(2);
      particle_kinds.push_back(2);
    };     
  };
    
  G4double m1 = dummy.getParticleMass(particle_kinds[0]);

  m1 *= m1;

  G4double m2 = dummy.getParticleMass(particle_kinds[1]);

  m2 *= m2;	 

  G4double a = 0.5 * (etot_scm * etot_scm - m1 - m2);

  G4double pmod = std::sqrt((a * a - m1 * m2) / (m1 + m2 + 2.0 * a));
  G4CascadeMomentum mom;
  std::pair<G4double, G4double> COS_SIN = randomCOS_SIN();
  G4double FI = randomPHI();
  G4double pt = pmod * COS_SIN.second;

  mom[1] = pt * std::cos(FI);
  mom[2] = pt * std::sin(FI);
  mom[3] = pmod * COS_SIN.first;

  G4CascadeMomentum mom1 = mom;

  for(G4int i = 1; i < 4; i++) mom1[i] *= -1.0;
  particles.push_back(G4InuclElementaryParticle(mom , particle_kinds[0]));
  particles.push_back(G4InuclElementaryParticle(mom1, particle_kinds[1]));

  return particles;
}

void G4ElementaryParticleCollider::initializeArrays()
{
  // Parameter array for momentum calculation of many body final states
  const G4double rmnData[14][10][2] = {
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

  // Copy to class scope
  for (G4int i = 0; i < 14; i++) {
    for (G4int j = 0; j < 10; j++) {
      for (G4int k = 0; k < 2; k++) rmn[i][j][k] = rmnData[i][j][k];
    }
  }

  // Parameter array for angular distribution calculation
  const G4double angData[4][4][13] = {
    {{ 2.7404, -30.853,  0.1026,-0.3829,  0.2499, 3.9025, 19.402, 
       0.1579, 0.3153,-17.953, 0.4217, 0.1499,  0.5369},
     {-9.6998,  106.24, -1.0542, 3.7587,  32.028,-91.126,-224.46, 
       2.9671,-7.4981, 109.72, 147.05, 2.8753, -13.216}, 
     { 10.400, -129.39,  11.389,-6.5144, -118.82, 323.73, 747.33,
       -5.5251, 43.295,-239.54,-653.35,-5.3078,  81.011}, 
     { 2.3882,  54.339, -16.638, 6.7740,  150.99,-400.48,-935.70, 
       6.8925,-76.460, 228.26, 915.07, 6.2233, -142.85}},

    {{-7.5137,  19.465, -0.4961, 103.81, -2.6994,-20.619,-44.180,
      -7.0218,-6.5373, 91.968,-3.5198,-5.9558, -10.550}, 
     { 44.096, -68.102,  11.800,-272.82, -460.45, 491.70, 471.94,
      -205.34, 193.07,-519.63,-260.19,-162.03,  296.29},
     {-74.379,  96.358, -90.857, 477.59,  1895.9,-1715.5,-1485.6, 
      569.51,-1018.1, 1126.6, 1225.0, 430.79, -1695.7}, 
     { 46.038, -56.827,  164.76,-512.22, -2519.0, 2114.3, 1805.5,
     -898.58, 1742.6,-1074.0,-1748.1,-625.48,  2893.5}},

    {{ 7.5479, -3.4831,  1.5437,-1788.2,  16.268, 33.004, 31.567, 
       134.96, 46.864,-132.70, 3.6373, 128.75,  69.621},
     {-39.274,  12.341, -33.769, 4305.2,  2138.4,-766.84,-301.76, 
       4872.2,-1303.0, 741.12, 155.92, 3140.2, -1924.5}, 
     { 64.835, -18.592,  251.92,-7931.4, -9126.2, 2700.3, 907.63,
       -14674., 6729.1,-1600.0,-752.01,-7918.9, 10620.0}, 
     { 41.609,  12.024, -450.71, 9347.1, 12431.0,-3352.5,-1077.3, 
       23924.,-11075., 1524.9, 1079.6, 10983., -17468.}},

    {{-1.8369,  0.1894, -1.2021, 7147.5, -29.654,-16.367,-6.8648,
      -821.16,-95.192, 58.598,-0.7804,-851.61, -138.65}, 
     { 8.6911, -0.6788,  0.2534,-3339.5, -3182.3, 373.94, 60.476,
      -32586., 2637.3,-318.74,-30.563,-18780.,  3928.1},
     {-13.060,  1.0665, -186.58,-4139.2,  13944.,-1320.2,-175.20,
       100980.,-12857., 677.51, 147.95, 44607., -20293.}, 
     { 7.1880, -0.7291,  332.54,-4436.4, -19342., 1642.3, 203.81,
       -165530.,20294.,-640.11,-212.50,-58790.,  32058.}}
  };

  // Copy to class scope
  for (G4int i = 0; i < 4; i++) {
    for (G4int j = 0; j < 4; j++) {
      for (G4int k = 0; k < 13; k++) ang[i][j][k] = angData[i][j][k];
    }
  }

  const G4double abnData[4][4][4] = {
    {{0.0856,  0.0716,  0.1729,  0.0376},  {5.0390,  3.0960,  7.1080,  1.4331},
     {-13.782, -11.125, -17.961, -3.1350},  {14.661,  18.130,  16.403,  6.4864}},
    {{0.0543,  0.0926, -0.1450,  0.2383}, {-9.2324, -3.2186, -13.032,  1.8253},
     {36.397,  20.273,  41.781,  1.7648}, {-42.962, -33.245, -40.799, -16.735}},
    {{-0.0511, -0.0515,  0.0454, -0.1541}, {4.6003,  0.8989,  8.3515, -1.5201},
     {-20.534, -7.5084, -30.260, -1.5692},  {27.731,  13.188,  32.882,  17.185}},
    {{0.0075,  0.0058, -0.0048,  0.0250}, {-0.6253, -0.0017, -1.4095,  0.3059},
     {2.9159,  0.7022,  5.3505,  0.3252}, {-4.1101, -1.4854, -6.0946, -3.5277}} 
  };

  // Copy to class scope
  for (G4int i = 0; i < 4; i++) {
    for (G4int j = 0; j < 4; j++) {
      for (G4int k = 0; k < 4; k++) abn[i][j][k] = abnData[i][j][k];
    }
  }

}
