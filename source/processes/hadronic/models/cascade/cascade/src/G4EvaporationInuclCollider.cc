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
#include "G4EvaporationInuclCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"
#include "G4ParticleLargerEkin.hh"
#include <algorithm>

typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;
	 
G4EvaporationInuclCollider::G4EvaporationInuclCollider()
  : verboseLevel(0) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::G4EvaporationInuclCollider" << G4endl;
  }
}

G4CollisionOutput G4EvaporationInuclCollider::collide(G4InuclParticle* /*bullet*/, G4InuclParticle* target) {

  verboseLevel = 0;
  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::evaporate" << G4endl;
  }

  G4CollisionOutput globalOutput;

  G4LorentzConvertor convertToTargetRestFrame;
  G4InuclNuclei* ntarget = dynamic_cast<G4InuclNuclei*>(target);
  convertToTargetRestFrame.setTarget(ntarget->getMomentum(), ntarget->getMass());

  G4double at = ntarget->getA();
  G4double zt = ntarget->getZ();
  G4double eEx = ntarget->getExitationEnergy();

  G4CascadeMomentum bmom;
  bmom[3] = convertToTargetRestFrame.getTRSMomentum();

  G4InuclNuclei targ(at, zt);
  G4CascadeMomentum tmom;
  targ.setExitationEnergy(eEx);
  targ.setMomentum(tmom);
  targ.setEnergy();

  targ.printParticle();

  G4CollisionOutput output;
  output = theEquilibriumEvaporator->collide(0, &targ);

  G4CollisionOutput TRFoutput;	   
  TRFoutput.addOutgoingParticles(output.getOutgoingParticles());  
  TRFoutput.addTargetFragments(output.getNucleiFragments());         

  if (verboseLevel > 3) {
    G4cout << " After EquilibriumEvaporator " << G4endl;
    output.printCollisionOutput();
  };
	 

  std::vector<G4InuclElementaryParticle> particles = TRFoutput.getOutgoingParticles();
  std::vector<G4InuclNuclei> nucleus = TRFoutput.getNucleiFragments();

  globalOutput.addOutgoingParticles(particles);
  globalOutput.addTargetFragments(nucleus);

  if (verboseLevel > 3) G4cout << "G4EvaporationInuclCollider::collide end" << G4endl;
 
	
  return globalOutput;
	
}
	
	     
G4bool G4EvaporationInuclCollider::inelasticInteractionPossible(G4InuclParticle* bullet,
						     G4InuclParticle* target, 
						     G4double ekin) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::inelasticInteractionPossible" << G4endl;
  }

  const G4double coeff = 0.001 * 1.2;
  const G4double one_third = 1.0 / 3.0;

  G4bool possible = true;
  G4double at;
  G4double zt;
  G4double ab;
  G4double zb;

  if (G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {
    at = nuclei_target->getA();
    zt = nuclei_target->getZ(); 
    if (G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {
      ab = nuclei_bullet->getA();
      zb = nuclei_bullet->getZ();     
    } else {
      G4InuclElementaryParticle* particle =
	dynamic_cast<G4InuclElementaryParticle*>(bullet);

      ab = 1;
      zb = particle->getCharge();
    }; 
  } else {
    if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {
      ab = nuclei_bullet->getA();
      zb = nuclei_bullet->getZ();     

      G4InuclElementaryParticle* particle =
	dynamic_cast<G4InuclElementaryParticle*>(target);

      at = 1;
      zt = particle->getCharge();    
    } else {

      return possible;
    };  
  }; 

  // VCOL used  for testing if elastic collision possible
  G4double VCOL = coeff * zt * zb / (std::pow(at, one_third) + std::pow(ab, one_third)); 

  // possible = VCOL < ekin; // NOTE: inelastic collision if not true
  possible = true; // we force elastic

  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::inelasticInteractionPossible" << G4endl;
    G4cout << " VCOL: " << VCOL << " ekin: " << ekin << " inelastic possible: " << possible << G4endl;
  }

  return possible;

}
	
G4InteractionCase G4EvaporationInuclCollider::bulletTargetSetter(G4InuclParticle* bullet,
						      G4InuclParticle* target) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::bulletTargetSetter" << G4endl;
  }

  G4InteractionCase interCase;

  if (G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {     
    if (G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) { // A + A         
      interCase.setInterCase(2);
      if (nuclei_target->getA() >= nuclei_bullet->getA()) {
	interCase.setBulletTarget(bullet, target);
      } else {
	interCase.setBulletTarget(target, bullet);
      }; 
    } else {
      interCase.setInterCase(1);
      interCase.setBulletTarget(bullet, target);
    }; 
  } else {
    G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet);
    if (nuclei_bullet) { 
      G4InuclElementaryParticle* part =
	dynamic_cast<G4InuclElementaryParticle*>(target);
      if (part) {
	interCase.setInterCase(1);
	interCase.setBulletTarget(target, bullet);
      };
    }; 
  };

  return interCase;
}       

G4bool G4EvaporationInuclCollider::explosion(G4InuclNuclei* target) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4EvaporationInuclCollider::explosion" << G4endl;
  }

  const G4double a_cut = 20.0;
  const G4double be_cut = 3.0;

  G4double a = target->getA();
  G4double z = target->getZ();
  G4double eexs = target->getExitationEnergy();
  G4bool explo = true;

  if (a > a_cut) {
    explo = false;
  } else {
    if (eexs < be_cut * bindingEnergy(a, z)) explo = false;
  };   

  return explo;
}
 
