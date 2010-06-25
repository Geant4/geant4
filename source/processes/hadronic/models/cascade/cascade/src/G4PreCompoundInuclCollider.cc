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
// $Id: G4PreCompoundInuclCollider.cc,v 1.10 2010-06-25 09:45:02 gunter Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Eliminate some unnecessary std::pow()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()"
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members, consolidate code

#include "G4PreCompoundInuclCollider.hh"
#include "G4BigBanger.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"
#include "G4NonEquilibriumEvaporator.hh"

	 
G4PreCompoundInuclCollider::G4PreCompoundInuclCollider()
  : G4VCascadeCollider("G4PreCompoundInuclCollider"),
    theElementaryParticleCollider(new G4ElementaryParticleCollider),
    theIntraNucleiCascader(new G4IntraNucleiCascader),
    theNonEquilibriumEvaporator(new G4NonEquilibriumEvaporator),
    theBigBanger(new G4BigBanger) {}

G4PreCompoundInuclCollider::~G4PreCompoundInuclCollider() {
  delete theElementaryParticleCollider;
  delete theIntraNucleiCascader;
  delete theNonEquilibriumEvaporator;
  delete theBigBanger;
}


void G4PreCompoundInuclCollider::collide(G4InuclParticle* bullet,
					 G4InuclParticle* target,
					 G4CollisionOutput& globalOutput) {
  if (verboseLevel > 3) {
    G4cout << " >>> G4PreCompoundInuclCollider::collide" << G4endl;
  }

  const G4int itry_max = 1000;
  		     
  if (useEPCollider(bullet,target)) {
    if (verboseLevel > 2) {
      bullet->printParticle();
      target->printParticle();
    }

    theElementaryParticleCollider->collide(bullet, target, globalOutput);
  } else { // needs to call all machinery    	
    G4LorentzConvertor convertToTargetRestFrame;

    interCase.set(bullet, target);
     
    if (interCase.valid()) { // ok
      G4InuclNuclei* ntarget =
	dynamic_cast<G4InuclNuclei*>(interCase.getTarget());

      convertToTargetRestFrame.setTarget(ntarget);
      G4int btype = 0;
      G4double ab = 0.0;
      G4double zb = 0.0;
      G4double at = ntarget->getA();
      G4double zt = ntarget->getZ();
       
      if (interCase.hadNucleus()) { // particle with nuclei
	G4InuclElementaryParticle* pbullet = 
	  dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());
         
	if (pbullet->isPhoton()) {
	  G4cout << " InuclCollider -> can not collide with photon " << G4endl;

	  globalOutput.trivialise(bullet, target);
	  return;
	} else {
	  convertToTargetRestFrame.setBullet(pbullet);   
	  btype = pbullet->type();
	}; 

      } else { // nuclei with nuclei
	G4InuclNuclei* nbullet = 
	  dynamic_cast<G4InuclNuclei*>(interCase.getBullet());

	convertToTargetRestFrame.setBullet(nbullet);   
	ab = nbullet->getA();
	zb = nbullet->getZ();
      };
       	
      G4double ekin = convertToTargetRestFrame.getKinEnergyInTheTRS();

      if (verboseLevel > 3) {
	G4cout << " ekin in trs " << ekin << G4endl;
      }

      if (inelasticInteractionPossible(bullet, target, ekin)) {
	convertToTargetRestFrame.toTheTargetRestFrame();

	if (verboseLevel > 3) {
	  G4cout << " degenerated? " << convertToTargetRestFrame.trivial() << G4endl;
	}

	G4LorentzVector bmom;
	bmom.setZ(convertToTargetRestFrame.getTRSMomentum());

	G4InuclNuclei ntarget(at, zt);		// Default is at rest

	theIntraNucleiCascader->setInteractionCase(interCase.code());
	 
	G4bool bad = true;
	G4int itry = 0;
	 
	G4CollisionOutput TRFoutput;
	G4CollisionOutput output;
	while (bad && itry < itry_max) {
	  itry++;

	  output.reset();	// Clear buffers for this attempt
	  TRFoutput.reset();

	  if (interCase.hadNucleus()) {
	    G4InuclElementaryParticle pbullet(bmom, btype);

	    theIntraNucleiCascader->collide(&pbullet, &ntarget, output);
	  } else {
	    G4InuclNuclei nbullet(bmom, ab, zb);
	    theIntraNucleiCascader->collide(&nbullet, &ntarget, output);
	  };   

	  if (verboseLevel > 3) {
	    G4cout << " After Cascade " << G4endl;
	    output.printCollisionOutput();
	  }
	  
	  // the rest, if any
	  TRFoutput.addOutgoingParticles(output.getOutgoingParticles());

	  if (output.numberOfNucleiFragments() == 1) { // there is smth. after
	    G4InuclNuclei cascad_rec_nuclei = output.getNucleiFragments()[0];
	    if (explosion(&cascad_rec_nuclei)) {
	      if (verboseLevel > 3) {
		G4cout << " big bang after cascade " << G4endl;
	      };

	      theBigBanger->collide(0,&cascad_rec_nuclei, TRFoutput);
	    } else {
	      output.reset();
	      theNonEquilibriumEvaporator->collide(0, &cascad_rec_nuclei, output);

	      if (verboseLevel > 3) {
		G4cout << " After NonEquilibriumEvaporator " << G4endl;
		output.printCollisionOutput();
	      };

	      TRFoutput.addOutgoingParticles(output.getOutgoingParticles());  
	      TRFoutput.addTargetFragments(output.getNucleiFragments());
	    };
	  };
	 
	  // convert to the LAB       
	  TRFoutput.boostToLabFrame(convertToTargetRestFrame);

	  globalOutput.addOutgoingParticles(TRFoutput.getOutgoingParticles());
	  globalOutput.addTargetFragments(TRFoutput.getNucleiFragments());
	  globalOutput.setOnShell(bullet, target);
	  if (globalOutput.acceptable()) return;

	  globalOutput.reset();		// Clear and try again
	};

	if (verboseLevel > 3) {
	  G4cout << " InuclCollider -> can not generate acceptable inter. after " 
		 << itry_max << " attempts " << G4endl;
	}
      } else {
	if (verboseLevel > 3) {
	  G4cout << " InuclCollider -> inelastic interaction is impossible " << G4endl
		 << " due to the coulomb barirer " << G4endl;
	}
      }
	
      globalOutput.trivialise(bullet, target);
      return;
    } else {
      if (verboseLevel > 3) {
	G4cout << " InuclCollider -> inter case " << interCase.code() << G4endl;
      };
    };       
  };

  return;
}
