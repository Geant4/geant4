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
// $Id: G4IntraNucleiCascader.cc,v 1.46 2010-07-03 00:07:55 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
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

#include "G4IntraNucleiCascader.hh"
#include "G4CascadParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4NucleiModel.hh"
#include "G4ParticleLargerEkin.hh"
#include "Randomize.hh"
#include <algorithm>

using namespace G4InuclSpecialFunctions;


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4IntraNucleiCascader::G4IntraNucleiCascader()
  : G4VCascadeCollider("G4IntraNucleiCascader") {}

G4IntraNucleiCascader::~G4IntraNucleiCascader() {}


void G4IntraNucleiCascader::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& output) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4IntraNucleiCascader::collide inter_case " << inter_case 
	   << G4endl;
  }

  G4NucleiModel model;
  model.setVerboseLevel(verboseLevel);
  theElementaryParticleCollider.setVerboseLevel(verboseLevel);

  const G4int itry_max = 1000;
  const G4int reflection_cut = 500;

  const G4double small_ekin = 0.001*MeV;	// Tolerance for round-off zero
  const G4double quasielast_cut = 1*MeV;	// To recover elastic scatters

  if (verboseLevel > 3) {
    bullet->printParticle();
    target->printParticle();
  }

  G4InuclNuclei* tnuclei = dynamic_cast<G4InuclNuclei*>(target);
  if (!tnuclei) {
    if (verboseLevel)
      G4cerr << " Target is not a nucleus.  Abandoning." << G4endl;
    return;
  }

  model.generateModel(tnuclei);

  G4double coulombBarrier = 0.00126*tnuclei->getZ()/
                                      (1.+G4cbrt(tnuclei->getA()));

  G4LorentzVector momentum_in = bullet->getMomentum() + target->getMomentum();

  // FIXME:  This assumes target at rest!  Should it be momentum_in.e()-m()?
  G4double ekin_in = bullet->getKineticEnergy();

  if (verboseLevel > 3) {
    model.printModel();
    G4cout << " intitial momentum  E " << momentum_in.e() << " Px "
	   << momentum_in.x() << " Py " << momentum_in.y() << " Pz "
	   << momentum_in.z() << G4endl;
  }

  // Bullet may be nucleus or simple particle
  G4InuclNuclei* bnuclei = dynamic_cast<G4InuclNuclei*>(bullet);
  G4InuclElementaryParticle* bparticle = 
                          dynamic_cast<G4InuclElementaryParticle*>(bullet);

  G4int itry = 0;
  while (itry < itry_max) {
    itry++;
    model.reset();

    std::vector<G4CascadParticle> cascad_particles;
    G4ExitonConfiguration theExitonConfiguration;
    std::vector<G4InuclElementaryParticle> output_particles;

    G4double afin = tnuclei->getA();	// Will deduct outgoing particles
    G4double zfin = tnuclei->getZ();	//    to determine recoil state
   
    if (inter_case == 1) { 		// particle with nuclei
      if (verboseLevel > 2)
	G4cout << " itry " << itry << " inter_case 1 " << G4endl;

      if (verboseLevel > 3)
	G4cout << " bparticle charge " << bparticle->getCharge()
	       << " baryon number " << bparticle->baryon() << G4endl;

      zfin += bparticle->getCharge();
      afin += bparticle->baryon();

      cascad_particles.push_back(model.initializeCascad(bparticle));
    } else {				// nuclei with nuclei
      if (verboseLevel > 2)
	G4cout << " itry " << itry << " inter_case " << inter_case << G4endl;

      G4double ab = bnuclei->getA();
      G4double zb = bnuclei->getZ();

      afin += ab;
      zfin += zb;

      G4NucleiModel::modelLists all_particles;    // Buffer to receive lists
      model.initializeCascad(bnuclei, tnuclei, all_particles);

      cascad_particles = all_particles.first;

      for (G4int ip = 0; ip < G4int(all_particles.second.size()); ip++) 
	output_particles.push_back(all_particles.second[ip]);

      if (cascad_particles.size() == 0) { // compound nuclei
	G4int ia = G4int(ab + 0.5);
	G4int iz = G4int(zb + 0.5);
	G4int i;

	for (i = 0; i < ia; i++) {
	  G4int knd = i < iz ? 1 : 2;
	  theExitonConfiguration.incrementQP(knd);
	};

	G4int ihn = G4int(2.0 * (ab - zb) * inuclRndm() + 0.5);
	G4int ihz = G4int(2.0 * zb * inuclRndm() + 0.5);

	for (i = 0; i < ihn; i++) theExitonConfiguration.incrementHoles(2);
	for (i = 0; i < ihz; i++) theExitonConfiguration.incrementHoles(1);
      };
    }; 

    std::vector<G4CascadParticle> new_cascad_particles;
    G4int iloop = 0;

    while (!cascad_particles.empty() && !model.empty()) {
      iloop++;

      if (verboseLevel > 2) {
	G4cout << " Number of cparticles " << cascad_particles.size() 
	       << " last one: " << G4endl;
	cascad_particles.back().print();
      }

      new_cascad_particles = model.generateParticleFate(cascad_particles.back(),
							&theElementaryParticleCollider);
      if (verboseLevel > 2) {
	G4cout << " New particles " << new_cascad_particles.size() << G4endl
	       << " Discarding last cparticle from list " << G4endl;
      }

      cascad_particles.pop_back();

      // handle the result of a new step

      if (new_cascad_particles.size() == 1) { // last particle goes without interaction
	if (model.stillInside(new_cascad_particles[0])) {
	  if (verboseLevel > 3) G4cout << " still inside " << G4endl;

	  if (new_cascad_particles[0].getNumberOfReflections() < reflection_cut &&
	      model.worthToPropagate(new_cascad_particles[0])) {
	    if (verboseLevel > 3) G4cout << " survives " << G4endl;
	    cascad_particles.push_back(new_cascad_particles[0]);
	  } else {
	    if (verboseLevel > 3) G4cout << " becomes an exiton " << G4endl;
	    theExitonConfiguration.incrementQP(new_cascad_particles[0].getParticle().type());
	  }
        } else { // particle about to leave nucleus - check for Coulomb barrier
	  if (verboseLevel > 3) G4cout << " possible escape " << G4endl;

          const G4InuclElementaryParticle& currentParticle =
	    new_cascad_particles[0].getParticle();

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
            if (KE > 0.0001) CBP = std::exp(-0.0181*0.5*tnuclei->getZ()*
                                            (1./KE - 1./coulombBarrier)*
                                         std::sqrt(mass*(coulombBarrier-KE)) );
            if (G4UniformRand() < CBP) {
	      if (verboseLevel > 3) {
		G4cout << " tunneled " << G4endl;
		currentParticle.printParticle();
	      }
	      output_particles.push_back(currentParticle);
            } else {
	      if (verboseLevel > 3) G4cout << " becomes an exciton " << G4endl;
              theExitonConfiguration.incrementQP(currentParticle.type());
            }
          } else {
	    if (verboseLevel > 3) {
	      G4cout << " Goes out " << G4endl;
	      currentParticle.printParticle();
	    }
	    output_particles.push_back(currentParticle);
          }
        } 
      } else { // interaction 
	if (verboseLevel > 3)
	  G4cout << " interacted, adding new to list " << G4endl;

	for (G4int i = 0; i < G4int(new_cascad_particles.size()); i++) 
	  cascad_particles.push_back(new_cascad_particles[i]);

	std::pair<G4int, G4int> holes = model.getTypesOfNucleonsInvolved();

	theExitonConfiguration.incrementHoles(holes.first);

	if (holes.second > 0)
	  theExitonConfiguration.incrementHoles(holes.second);
      }		// if (new_cascad_particles ...
    }		// while cascade-list and model

    // Cascade is finished. Check if it's OK.

    if (verboseLevel > 3) {
      G4cout << " Cascade finished  " << G4endl
	     << " output_particles  " << output_particles.size() <<  G4endl;
    }

    G4LorentzVector momentum_out;
    particleIterator ipart;

    for (ipart = output_particles.begin(); ipart != output_particles.end(); ipart++) {
      if (verboseLevel > 3) {
	ipart->printParticle();
	G4cout << "  charge " << ipart->getCharge() << " baryon number "
	       << ipart->baryon() << G4endl;
      }

      momentum_out += ipart->getMomentum();

      zfin -= ipart->getCharge();
      afin -= ipart->baryon();
    };

    if (afin<0. || zfin<0. || afin<zfin) {  // Sanity check before proceeding
      G4cerr << " >>> G4IntraNucleiCascader ERROR:  Recoil nucleus is not"
	     << " physical! A=" << afin << " Z=" << zfin << G4endl;
      continue;				// Discard event and try again
    }

    G4LorentzVector presid = momentum_in - momentum_out;

    if (verboseLevel > 3) {
      G4cout << "  afin " << afin << " zfin " << zfin <<  G4endl;

      G4cout << " momentum_in:    px " << momentum_in.px()
	     << " py " << momentum_in.py() << " pz " << momentum_in.pz()
	     << " E " << momentum_in.e() << G4endl;
      G4cout << " momentum_out:   px " << momentum_out.px()
	     << " py " << momentum_out.py() << " pz " << momentum_out.pz()
	     << " E " << momentum_out.e() << G4endl;
      G4cout << " resid (in-out): px " << presid.px()
	     << " py " << presid.py() << " pz " << presid.pz()
	     << " E " << presid.e() << G4endl;
    }

    if (afin > 1.0) {
      // Excitation is difference between outgoing and ground-state mass
      G4double mass = G4InuclNuclei::getNucleiMass(afin, zfin);
      G4double mres = presid.m();

      G4double Eex = (mres-mass)*GeV;		// Excitation is in MeV
      if (std::abs(Eex) < small_ekin) Eex = 0.;	// Absorb round-off

      // Quasi-elastic scattering 
      if (output_particles.size() == 1 && std::abs(Eex) < quasielast_cut) {
	if (verboseLevel > 3)
	  G4cout << " quasi-elastic scatter with " << Eex << " MeV recoil"
		 << G4endl;
	Eex = 0.;
      }

      if (Eex < 0.0) {		// Unphysical recoil; reject cascade
	if (verboseLevel > 2) G4cout << " unphysical recoil " << Eex << G4endl;
	continue;
      }

      if (verboseLevel > 3) {
	G4cout << " fragment mass " << mass << " recoil " << mres
	       << "  Eex  " << Eex  << " MeV" << G4endl;
      }

      if (!goodCase(afin, zfin, Eex, ekin_in)) { 	// unphysical recoil
	if (verboseLevel > 2) G4cout << " nuclear recoil failed." << G4endl;
	continue;
      }

      if (verboseLevel > 2)
	G4cout << " adding recoil nucleus/fragment to output list" << G4endl;

      G4InuclNuclei outgoing_nuclei(presid, afin, zfin, Eex);
      outgoing_nuclei.setModel(4);
      outgoing_nuclei.setExitonConfiguration(theExitonConfiguration);
      if (verboseLevel > 3) outgoing_nuclei.printParticle();
      
      output.addTargetFragment(outgoing_nuclei);

      std::sort(output_particles.begin(), output_particles.end(), G4ParticleLargerEkin());
      output.addOutgoingParticles(output_particles);
      
      return;
    } else { 	// special case, when one has no nuclei after the cascad
      if (afin == 1.0) { // recoiling nucleon
	G4int last_type = (zfin == 1.) ? 1 : 2;

	G4double mass = G4InuclElementaryParticle::getParticleMass(last_type);
	G4double mres = presid.m();

	if (mres-mass < -small_ekin) {		// Below round-off tolerance
	  if (verboseLevel > 2) 
	    G4cout << " unphysical recoil nucleon" << G4endl;
	  continue;
	}

	if (verboseLevel > 3)
	  G4cout << " adding recoiling nucleon to output list" << G4endl;

	G4InuclElementaryParticle last_particle(presid, last_type, 4);
	if (verboseLevel > 3) last_particle.printParticle();
	
	output_particles.push_back(last_particle);
      }; 

      std::sort(output_particles.begin(), output_particles.end(), G4ParticleLargerEkin());
      output.addOutgoingParticles(output_particles);

      return;
    }		// if (afin > 1)
  }		// while (itry < itry_max)

  if (verboseLevel > 3) {
    G4cout << " IntraNucleiCascader-> no inelastic interaction after "
	   << itry_max << " attempts " << G4endl;
  }

  output.trivialise(bullet, target);

  return;
}

G4bool G4IntraNucleiCascader::goodCase(G4double a, 
				       G4double z, 
				       G4double eexs, 
				       G4double ein) const {
  if (verboseLevel > 1) {
    G4cout << " >>> G4IntraNucleiCascader::goodCase" << G4endl;
  }

  const G4double eexs_cut = 0.0001*MeV;
  const G4double reason_cut = 7.0*MeV;
  const G4double ediv_cut = 5.0*MeV;

  G4bool good = (eexs <= eexs_cut);	// Always allow for unexcited recoil

  if (eexs > eexs_cut) {
    G4double eexs_max0z = ein*GeV / ediv_cut;   // ein is GeV, eexs is MeV
    G4double dm = bindingEnergy(a,z);
    G4double eexs_max = eexs_max0z > reason_cut*dm ? eexs_max0z : reason_cut * dm;

    if (verboseLevel > 3) {
      G4cout << " eexs " << eexs << " max " << eexs_max << " dm " << dm << G4endl;
    }

    good = (eexs < eexs_max);
  };

  return good; 
}
