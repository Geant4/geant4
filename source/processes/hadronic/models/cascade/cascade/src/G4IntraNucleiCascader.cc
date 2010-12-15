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
// $Id: G4IntraNucleiCascader.cc,v 1.69 2010-12-15 07:41:09 gunter Exp $
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

#include "G4IntraNucleiCascader.hh"
#include "G4CascadParticle.hh"
#include "G4CascadeRecoilMaker.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4CollisionOutput.hh"
#include "G4DecayTable.hh"
#include "G4DecayProducts.hh"
#include "G4HadTmpUtil.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclSpecialFunctions.hh"
#include "G4LorentzConvertor.hh"
#include "G4NucleiModel.hh"
#include "G4ParticleLargerEkin.hh"
#include "Randomize.hh"
#include <algorithm>

using namespace G4InuclSpecialFunctions;


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4IntraNucleiCascader::G4IntraNucleiCascader()
  : G4CascadeColliderBase("G4IntraNucleiCascader"),
    model(new G4NucleiModel),
    theElementaryParticleCollider(new G4ElementaryParticleCollider),
    theRecoilMaker(new G4CascadeRecoilMaker) {}

G4IntraNucleiCascader::~G4IntraNucleiCascader() {
  delete model;
  delete theElementaryParticleCollider;
  delete theRecoilMaker;
}


void G4IntraNucleiCascader::collide(G4InuclParticle* bullet,
				    G4InuclParticle* target,
				    G4CollisionOutput& globalOutput) {
  if (verboseLevel) G4cout << " >>> G4IntraNucleiCascader::collide " << G4endl;

  const G4int itry_max = 1000;
  const G4int reflection_cut = 500;

  const G4double small_ekin = 0.001*MeV;	// Tolerance for round-off zero
  const G4double quasielast_cut = 1*MeV;	// To recover elastic scatters

  // Configure processing modules
  model->setVerboseLevel(verboseLevel);
  theElementaryParticleCollider->setVerboseLevel(verboseLevel);
  theRecoilMaker->setVerboseLevel(verboseLevel);
  theRecoilMaker->setTolerance(small_ekin);

  // Energy/momentum conservation usually requires a recoiling nuclear fragment
  // This cut will be increased on each "itry" if momentum could not balance.
  G4double minimum_recoil_A = 0.;		// Nuclear fragment required

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

  interCase.set(bullet,target);		// Classify collision type

  model->generateModel(tnuclei);

  G4double coulombBarrier = 0.00126*tnuclei->getZ()/
                                      (1.+G4cbrt(tnuclei->getA()));

  G4LorentzVector momentum_in = bullet->getMomentum() + target->getMomentum();

  if (verboseLevel > 3) {
    model->printModel();
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
    if (verboseLevel > 2) {
      G4cout << " itry " << itry << " inter_case " << interCase.code()
	     << G4endl;
    }

    model->reset();    			// Start new cascade process
    output.reset();
    cascad_particles.clear();
    output_particles.clear();
    theExitonConfiguration.clear();

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

      output_particles.insert(output_particles.end(),
			      all_particles.second.begin(),
			      all_particles.second.end());

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

    new_cascad_particles.clear();
    G4int iloop = 0;

    while (!cascad_particles.empty() && !model->empty()) {
      iloop++;

      if (verboseLevel > 2) {
	G4cout << " Iteration " << iloop << ": Number of cparticles "
	       << cascad_particles.size() << " last one: " << G4endl;
	cascad_particles.back().print();
      }

      new_cascad_particles = model->generateParticleFate(cascad_particles.back(),
							theElementaryParticleCollider);
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
            if (KE > 0.0001) CBP = std::exp(-0.0181*0.5*tnuclei->getZ()*
                                            (1./KE - 1./coulombBarrier)*
                                         std::sqrt(mass*(coulombBarrier-KE)) );

            if (G4UniformRand() < CBP) {
	      if (verboseLevel > 3) {
		G4cout << " tunneled " << G4endl;
		currentParticle.printParticle();
	      }
	      // Tunnelling through barrier leaves KE unchanged
	      output_particles.push_back(currentParticle);
            } else {
	      processTrappedParticle(currentCParticle);
            }
          } else {
	    if (verboseLevel > 3) G4cout << " Goes out " << G4endl;

	    output_particles.push_back(currentParticle);

	    /*****
	    // Adjust kinetic energy by height of potential (+ve or -ve)
	    G4double newKE = KE - Q*coulombBarrier;
	    output_particles.back().setKineticEnergy(newKE);
	    *****/

	    if (verboseLevel > 3) output_particles.back().printParticle();
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
      theRecoilMaker->collide(bullet,target,output_particles,cascad_particles);

      G4double aresid = theRecoilMaker->getRecoilA();
      if (verboseLevel > 2) {
	G4cout << " cparticles remaining " << cascad_particles.size()
	       << " nucleus (model) has "
	       << model->getNumberOfNeutrons() << " n, "
	       << model->getNumberOfProtons() << " p "
	       << " residual fragment A " << aresid << G4endl;
      }

      if (aresid <= minimum_recoil_A) break;	// Must have minimum fragment
    }		// while cascade-list and model

    // Add left-over cascade particles to output
    for (G4int i = 0; i < G4int(cascad_particles.size()); i++)
      output_particles.push_back(cascad_particles[i].getParticle());
 
    // Cascade is finished. Check if it's OK.
    if (verboseLevel > 3) {
      G4cout << " Cascade finished  " << G4endl
	     << " output_particles  " << output_particles.size() <<  G4endl;

      particleIterator ipart = output_particles.begin();
      for (; ipart != output_particles.end(); ipart++) {
	ipart->printParticle();
	G4cout << "  charge " << ipart->getCharge() << " baryon number "
	       << ipart->baryon() << G4endl;
      }
    }

    // Use last created recoil fragment instead of re-constructing
    G4int afin = theRecoilMaker->getRecoilA();
    G4int zfin = theRecoilMaker->getRecoilZ();

    // Sanity check before proceeding
    if (!theRecoilMaker->goodFragment() && !theRecoilMaker->wholeEvent()) {
      if (verboseLevel > 1)
	G4cerr << " Recoil nucleus is not physical: A=" << afin << " Z="
	       << zfin << G4endl;
      continue;				// Discard event and try again
    }

    const G4LorentzVector& presid = theRecoilMaker->getRecoilMomentum();

    if (verboseLevel > 1) {
      G4cout << "  afin " << afin << " zfin " << zfin <<  G4endl;
    }

    if (afin == 0) break;		// Whole event fragmented, exit

    if (afin == 1) {			// Add bare nucleon to particle list
      G4int last_type = (zfin==1) ? 1 : 2;	// proton=1, neutron=2

      G4double mass = G4InuclElementaryParticle::getParticleMass(last_type);
      G4double mres = presid.m();

      // Check for sensible kinematics
      if (mres-mass < -small_ekin) {		// Insufficient recoil energy
	if (verboseLevel > 2) G4cerr << " unphysical recoil nucleon" << G4endl;
	continue;
      }

      if (mres-mass > small_ekin) {		// Too much extra energy
	if (verboseLevel > 2)
	  G4cerr << " extra energy with recoil nucleon" << G4endl;

	// FIXME:  For now, we add the nucleon as unbalanced, and let
	//	   "SetOnShell" fudge things.  This should be abandoned.
      }

      G4InuclElementaryParticle last_particle(presid, last_type, 4);

      if (verboseLevel > 3) {
	G4cout << " adding recoiling nucleon to output list" << G4endl;
	last_particle.printParticle();
      }

      output_particles.push_back(last_particle);
    }

    // Process recoil fragment for consistency, exit or reject
    if (output_particles.size() == 1) {
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
	continue;
      }
      output.addRecoilFragment(*recoilFrag);

      // TEMPORARY:  Add both frag and nuclei, for code validation
      G4InuclNuclei* recoilNucl = theRecoilMaker->makeRecoilNuclei(4);
      if (!recoilFrag) {
	G4cerr << "Got null pointer for recoil nucleus!" << G4endl;
	continue;
      }
      output.addOutgoingNucleus(*recoilNucl);
      
      if (verboseLevel > 2)
	G4cout << " adding recoil nucleus/fragment to output list" << G4endl;
    }

    // Put final-state particle in "leading order" for return
    std::sort(output_particles.begin(), output_particles.end(), G4ParticleLargerEkin());
    output.addOutgoingParticles(output_particles);

    // Adjust final state without fragment to balance momentum and energy
    if (afin <= 1) {
      output.setVerboseLevel(verboseLevel);
      output.setOnShell(bullet, target);
      output.setVerboseLevel(0);

      if (output.acceptable()) break;
    } else if (theRecoilMaker->goodNucleus()) break;

    // Cascade not physically reasonable
    if (afin <= minimum_recoil_A && minimum_recoil_A < tnuclei->getA()) {
      ++minimum_recoil_A;
      if (verboseLevel > 3) {
	G4cout << " minimum recoil fragment increased to A " << minimum_recoil_A
	       << G4endl;
      }
    }
  }	// while (itry < itry_max)

  // Cascade completed, for good or ill
  if (itry == itry_max) {
    if (verboseLevel > 3) {
      G4cout << " IntraNucleiCascader-> no inelastic interaction after "
	     << itry_max << " attempts " << G4endl;
    }

    output.trivialise(bullet, target);
  } else if (verboseLevel) {
    G4cout << " IntraNucleiCascader output after trials " << itry << G4endl;
  }

  // Copy final generated cascade to output buffer for return
  globalOutput.add(output);
  return;
}


// Convert particles which cannot escape into excitons (or eject/decay them)

void G4IntraNucleiCascader::
processTrappedParticle(const G4CascadParticle& trapped) {
  const G4InuclElementaryParticle& trappedP = trapped.getParticle();

  G4int xtype = trappedP.type();
  if (verboseLevel > 3) G4cout << " exciton of type " << xtype << G4endl;
  
  if (trappedP.nucleon()) {	// normal exciton (proton or neutron)
    theExitonConfiguration.incrementQP(xtype);
    return;
  }

  if (trappedP.hyperon()) {	// Not nucleon, so must be hyperon
    decayTrappedParticle(trapped);
    return;
  }

  // non-standard exciton; release it
  // FIXME: this is a meson, so need to absorb it
  if (verboseLevel > 3) {
    G4cout << " non-standard should be absorbed, now released" << G4endl;
    trapped.print();
  }
  
  output_particles.push_back(trappedP);
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

    output_particles.push_back(trappedP);
    return;
  }

  // Get secondaries from decay in particle's rest frame
  G4DecayProducts* daughters = unstable->SelectADecayChannel()->DecayIt();
  if (!daughters) {			// No final state; cannot decay!
    if (verboseLevel > 3)
      G4cerr << " no daughters from trapped particle decay" << G4endl;

    output_particles.push_back(trappedP);
    return;
  }

  if (verboseLevel > 3) daughters->DumpInfo();

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

    G4InuclElementaryParticle idaugEP(*idaug, 4);

    // Only hadronic secondaries can be propagated; photons escape
    if (idaugEP.isPhoton()) output_particles.push_back(idaugEP);
    else {
      G4CascadParticle idaugCP(idaugEP, decayPos, zone, 0., gen);
      cascad_particles.push_back(idaugCP);
    }
  }
}
