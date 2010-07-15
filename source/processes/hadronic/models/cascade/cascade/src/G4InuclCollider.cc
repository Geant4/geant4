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
// $Id: G4InuclCollider.cc,v 1.40 2010-07-15 19:34:09 mkelsey Exp $
// Geant4 tag: $Name: not supported by cvs2svn $
//
// 20100114  M. Kelsey -- Remove G4CascadeMomentum, use G4LorentzVector directly
// 20100309  M. Kelsey -- Eliminate some unnecessary std::pow()
// 20100413  M. Kelsey -- Pass G4CollisionOutput by ref to ::collide()
// 20100418  M. Kelsey -- Move lab-frame transformation code to G4CollisonOutput
// 20100429  M. Kelsey -- Change "photon()" to "isPhoton()"
// 20100517  M. Kelsey -- Inherit from common base class, make other colliders
//		simple data members, consolidate code
// 20100620  M. Kelsey -- Reorganize top level if-blocks to reduce nesting,
//		use new four-vector conservation check.
// 20100701  M. Kelsey -- Bug fix energy-conservation after equilibrium evap,
//		pass verbosity through to G4CollisionOutput
// 20100714  M. Kelsey -- Move conservation checking to base class, report
//		number of iterations at end

#include "G4InuclCollider.hh"
#include "G4BigBanger.hh"
#include "G4CollisionOutput.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"
#include "G4NonEquilibriumEvaporator.hh"


G4InuclCollider::G4InuclCollider()
  : G4CascadeColliderBase("G4InuclCollider") {}

G4InuclCollider::~G4InuclCollider() {}


void G4InuclCollider::collide(G4InuclParticle* bullet, G4InuclParticle* target,
			      G4CollisionOutput& globalOutput) {
  if (verboseLevel > 1) {
    G4cout << " >>> G4InuclCollider::collide" << G4endl;
  }

  // Initialize colliders verbosity
  theElementaryParticleCollider.setVerboseLevel(verboseLevel);
  theIntraNucleiCascader.setVerboseLevel(verboseLevel);
  theNonEquilibriumEvaporator.setVerboseLevel(verboseLevel);
  theEquilibriumEvaporator.setVerboseLevel(verboseLevel);
  theBigBanger.setVerboseLevel(verboseLevel);

  TRFoutput.setVerboseLevel(verboseLevel);
  output.setVerboseLevel(verboseLevel);

  const G4int itry_max = 1000;

  // Particle-on-particle collision; no nucleus involved
  if (useEPCollider(bullet,target)) {
    if (verboseLevel > 2)
      G4cout << " InuclCollider -> particle on particle collision" << G4endl;
 
    theElementaryParticleCollider.collide(bullet, target, globalOutput);
    return;
  }
  
  interCase.set(bullet,target);		// Classify collision type
  if (verboseLevel > 2) {
    G4cout << " InuclCollider -> inter case " << interCase.code() << G4endl;
  }

  if (!interCase.valid()) {
    if (verboseLevel > 1)
      G4cerr << " InuclCollider -> no collision possible " << G4endl;

    globalOutput.trivialise(bullet, target);
    return;
  }

  // Target must be a nucleus
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(interCase.getTarget());
    
  G4LorentzConvertor convertToTargetRestFrame;
  convertToTargetRestFrame.setTarget(ntarget);
  G4int btype = 0;
  G4double ab = 0.0;
  G4double zb = 0.0;
  
  if (interCase.hadNucleus()) { 	// particle with nuclei
    G4InuclElementaryParticle* pbullet = 
      dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());
    
    if (pbullet->isPhoton()) {
      G4cerr << " InuclCollider -> can not collide with photon " << G4endl;
      globalOutput.trivialise(bullet, target);
      return;
    } else {
      convertToTargetRestFrame.setBullet(pbullet);   
      btype = pbullet->type();
    } 
  } else { 				// nuclei with nuclei
    G4InuclNuclei* nbullet = 
      dynamic_cast<G4InuclNuclei*>(interCase.getBullet());
    
    convertToTargetRestFrame.setBullet(nbullet);   
    ab = nbullet->getA();
    zb = nbullet->getZ();
  }
  
  G4double ekin = convertToTargetRestFrame.getKinEnergyInTheTRS();
  
  if (verboseLevel > 3) G4cout << " ekin in trs " << ekin << G4endl;

  if (!inelasticInteractionPossible(bullet, target, ekin)) {
    if (verboseLevel > 3)
      G4cout << " InuclCollider -> inelastic interaction is impossible\n"
	     << " due to the coulomb barirer " << G4endl;

    globalOutput.trivialise(bullet, target);
    return;
  }

  // Generate interaction secondaries in rest frame of target nucleus
  convertToTargetRestFrame.toTheTargetRestFrame();
  if (verboseLevel > 3) {
    G4cout << " degenerated? " << convertToTargetRestFrame.trivial()
	   << G4endl;
  }
  
  G4LorentzVector bmom;			// Bullet is along local Z
  bmom.setZ(convertToTargetRestFrame.getTRSMomentum());
  
  theIntraNucleiCascader.setInteractionCase(interCase.code());

  G4int itry = 0;
  while (itry < itry_max) {
    itry++;
    if (verboseLevel > 2)
      G4cout << " IntraNucleiCascader itry " << itry << G4endl;

    globalOutput.reset();		// Clear buffers for this attempt
    output.reset();	
    TRFoutput.reset();
    
    if (interCase.hadNucleus()) {
      G4InuclElementaryParticle pbullet(bmom, btype);
      theIntraNucleiCascader.collide(&pbullet, ntarget, output);
    } else {
      G4InuclNuclei nbullet(bmom, ab, zb);
      theIntraNucleiCascader.collide(&nbullet, ntarget, output);
    }   
    
    if (verboseLevel > 1) {
      G4cout << " After Cascade " << G4endl;
    }

    // the rest, if any
    // FIXME:  The code below still does too much copying!
    TRFoutput.addOutgoingParticles(output.getOutgoingParticles());
    
    if (output.numberOfNucleiFragments() == 1) { // there is smth. after
      G4InuclNuclei cascad_rec_nuclei = output.getNucleiFragments()[0];
      if (explosion(&cascad_rec_nuclei)) {
	if (verboseLevel > 1) {
	  G4cout << " big bang after cascade " << G4endl;
	}

	// Add result of explosion diretly to output
	theBigBanger.collide(0,&cascad_rec_nuclei, TRFoutput);
      } else {
	output.reset();			// Buffer for debugging evaporators
	theNonEquilibriumEvaporator.collide(0, &cascad_rec_nuclei, output);
	
	if (verboseLevel > 1) {
	  G4cout << " After NonEquilibriumEvaporator " << G4endl;
	}

	TRFoutput.addOutgoingParticles(output.getOutgoingParticles());

	// Use nuclear fragment left from non-equilibrium for next step
	G4InuclNuclei exiton_rec_nuclei = output.getNucleiFragments()[0];
	
	output.reset();
	theEquilibriumEvaporator.collide(0, &exiton_rec_nuclei, output);
	
	if (verboseLevel > 1) {
	  G4cout << " After EquilibriumEvaporator " << G4endl;
	}

	TRFoutput.addOutgoingParticles(output.getOutgoingParticles());  
	TRFoutput.addTargetFragments(output.getNucleiFragments());
      }
    }

    if (verboseLevel > 2)
      G4cout << " itry " << itry << " finished, moving to lab frame" << G4endl;

    // convert to the LAB
    TRFoutput.boostToLabFrame(convertToTargetRestFrame);
    globalOutput.addOutgoingParticles(TRFoutput.getOutgoingParticles());
    globalOutput.addTargetFragments(TRFoutput.getNucleiFragments());

    // Check energy conservation before mucking everything up
    balance.setOwner("Before-setOnShell");
    validateOutput(bullet, target, globalOutput);
    balance.setOwner(theName);

    // Adjust final state particles to balance momentum and energy
    globalOutput.setOnShell(bullet, target);
    if (globalOutput.acceptable()) {
      if (verboseLevel) 
	G4cout << " InuclCollider output after trials " << itry << G4endl;

      return;
    }
  }
  
  if (verboseLevel) {
    G4cout << " InuclCollider -> can not generate acceptable inter. after " 
	   << itry_max << " attempts " << G4endl;
  }
  
  globalOutput.trivialise(bullet, target);
  return;
}
