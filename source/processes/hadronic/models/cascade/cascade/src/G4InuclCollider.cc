#include "G4InuclCollider.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4LorentzConvertor.hh"
#include "G4ParticleLargerEkin.hh"
#include "algorithm"

typedef vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef vector<G4InuclNuclei>::iterator nucleiIterator;
	 
G4InuclCollider::G4InuclCollider()
  : verboseLevel(2) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclCollider::G4InuclCollider" << G4endl;
  }
}

G4CollisionOutput G4InuclCollider::collide(G4InuclParticle* bullet,
					   G4InuclParticle* target) {

  verboseLevel = 2;
  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclCollider::collide" << G4endl;
  }

  const G4int itry_max = 1000;
  		     
  G4CollisionOutput globalOutput;
  G4InuclElementaryParticle* particle1 =
    dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4InuclElementaryParticle* particle2 =
    dynamic_cast<G4InuclElementaryParticle*>(target);
  
  if(particle1 && particle2) { // particle + particle 
    if (verboseLevel > 1) {
      particle1->printParticle();
      particle2->printParticle();
    }
    globalOutput = theElementaryParticleCollider->collide(bullet, target);
  }
  else { // needs to call all machinery    	

    G4LorentzConvertor convertToTargetRestFrame;
    G4InteractionCase interCase = bulletTargetSetter(bullet, target);
    G4int intcase = interCase.getInterCase();
     
    if(intcase > 0) { // ok

      G4InuclNuclei* ntarget =
	dynamic_cast<G4InuclNuclei*>(interCase.getTarget());

      convertToTargetRestFrame.setTarget(ntarget->getMomentum(),
					 ntarget->getMass());

      G4int btype = 0;
      G4double ab = 0.0;
      G4double zb = 0.0;
      G4double at = ntarget->getA();
      G4double zt = ntarget->getZ();
       
      if(intcase == 1) { // particle with nuclei

	G4InuclElementaryParticle* pbullet = 
	  dynamic_cast<G4InuclElementaryParticle*>(interCase.getBullet());
         
	if(pbullet->photon()) {

	  G4cout << " InuclCollider -> can not collide with photon " << G4endl;

	  globalOutput.trivialise(bullet, target);

	  return globalOutput;
	}
	else {
	  convertToTargetRestFrame.setBullet(pbullet->getMomentum(),
					     pbullet->getMass());   
	  btype = pbullet->type();
	}; 
      }
      else { // nuclei with nuclei

	G4InuclNuclei* nbullet = 
	  dynamic_cast<G4InuclNuclei*>(interCase.getBullet());

	convertToTargetRestFrame.setBullet(nbullet->getMomentum(),
					   nbullet->getMass());   
	ab = nbullet->getA();
	zb = nbullet->getZ();
      };
       	
      G4double ekin = convertToTargetRestFrame.getKinEnergyInTheTRS();

      if (verboseLevel > 1) {
	G4cout << " ekin in trs " << ekin << G4endl;
      }

      if(inelasticInteractionPossible(bullet, target, ekin)) {
	convertToTargetRestFrame.toTheTargetRestFrame();

	if (verboseLevel > 1) {
	  G4cout << " degenerated? " << convertToTargetRestFrame.trivial() << G4endl;
	}

	vector<G4double> bmom(4, 0.0);

	bmom[3] = convertToTargetRestFrame.getTRSMomentum();

	G4InuclNuclei ntarget(at, zt);
	vector<G4double> tmom(4, 0.0);

	ntarget.setMomentum(tmom);
	ntarget.setEnergy();
	theIntraNucleiCascader->setInteractionCase(intcase);
	 
	G4bool bad = true;
	G4int itry = 0;
	 
	while (bad && itry < itry_max) {

	  G4CollisionOutput TRFoutput;
	  G4CollisionOutput output;

	  itry++;
	  if(intcase == 1) {

	    G4InuclElementaryParticle pbullet(bmom, btype);

	    output = theIntraNucleiCascader->collide(&pbullet, &ntarget);
	  }
	  else {

	    G4InuclNuclei nbullet(ab, zb);

	    nbullet.setMomentum(bmom);
	    nbullet.setEnergy();
	    output = theIntraNucleiCascader->collide(&nbullet, &ntarget);
	  };   

	  if (verboseLevel > 0) {
	    G4cout << " After Cascade " << G4endl;

	    output.printCollisionOutput();

	    G4cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << 
	      G4endl;
	  }
	  
	  // the rest, if any
	  TRFoutput.addOutgoingParticles(output.getOutgoingParticles());

	  if(output.numberOfNucleiFragments() == 1) { // there is smth. after
	    
	    G4cout << "::::::::::::::" << G4endl;
	    
	    G4InuclNuclei cascad_rec_nuclei = output.getNucleiFragments()[0];

	    if(explosion(&cascad_rec_nuclei)) {

	      G4cout << " big bang after cascade " << G4endl;

	      output = theBigBanger->collide(0,&cascad_rec_nuclei);
	      TRFoutput.addOutgoingParticles(output.getOutgoingParticles());
	    }
	    else {
	      output = theNonEquilibriumEvaporator->collide(0, &cascad_rec_nuclei);

	      if (verboseLevel > 1) {
		G4cout << " After NonEquilibriumEvaporator " << G4endl;

		output.printCollisionOutput();

		G4cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << G4endl;
	      }

	      TRFoutput.addOutgoingParticles(output.getOutgoingParticles());

	      G4InuclNuclei exiton_rec_nuclei = output.getNucleiFragments()[0];

	      output = theEquilibriumEvaporator->collide(0, &exiton_rec_nuclei);

	      if (verboseLevel > 1) {
		G4cout << " After EquilibriumEvaporator " << G4endl;

		output.printCollisionOutput();

		G4cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << G4endl;
	      }

	      TRFoutput.addOutgoingParticles(output.getOutgoingParticles());  
	      TRFoutput.addTargetFragments(output.getNucleiFragments());         
	    };
	  };

	 
	  // convert to the LAB       
	  G4bool withReflection = convertToTargetRestFrame.reflectionNeeded();       
	  vector<G4InuclElementaryParticle> particles = 
	    TRFoutput.getOutgoingParticles();

	  if(!particles.empty()) { 
	    particleIterator ipart;
	    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {

	      vector<G4double> mom = ipart->getMomentum();

	      if(withReflection) mom[3] = -mom[3];
	      mom = convertToTargetRestFrame.rotate(mom);
	      ipart->setMomentum(mom); 
	      mom = convertToTargetRestFrame.backToTheLab(ipart->getMomentum());
	      ipart->setMomentum(mom); 
	    };
	    sort(particles.begin(), particles.end(), G4ParticleLargerEkin());
	  };
           
	  vector<G4InuclNuclei> nucleus = TRFoutput.getNucleiFragments();

	  if(!nucleus.empty()) { 
	    nucleiIterator inuc;
	    for(inuc = nucleus.begin(); inuc != nucleus.end(); inuc++) {
	      vector<G4double> mom = inuc->getMomentum(); 
	      if(withReflection) mom[3] = -mom[3];
	      mom = convertToTargetRestFrame.rotate(mom);
	      inuc->setMomentum(mom);
	      inuc->setEnergy(); 
	      mom = convertToTargetRestFrame.backToTheLab(inuc->getMomentum());
	      inuc->setMomentum(mom);
	      inuc->setEnergy(); 
	    };
	  };
	  globalOutput.addOutgoingParticles(particles);
	  globalOutput.addTargetFragments(nucleus);
	  globalOutput.setOnShell(bullet, target);
	  if(globalOutput.acceptable()) {

	    return globalOutput;
	  }
	  else {
	    globalOutput.reset();
	  }; 
	};

	if (verboseLevel > 1) {
	  G4cout << " InuclCollider -> can not generate acceptable inter. after " 
		 << itry_max << " attempts " << G4endl;
	}

	globalOutput.trivialise(bullet, target);

	return globalOutput;        
      }
      else {

	if (verboseLevel > 1) {
	  G4cout << " InuclCollider -> inelastic interaction is impossible " << G4endl
		 << " due to the coulomb barirer " << G4endl;
	}

	globalOutput.trivialise(bullet, target);

	return globalOutput;
      };	
    }
    else {

      G4cout << " InuclCollider -> inter case " << intcase << G4endl;

    };       
  };

  return globalOutput;
}
		     
G4bool G4InuclCollider::inelasticInteractionPossible(G4InuclParticle* bullet,
						     G4InuclParticle* target, 
						     G4double ekin) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclCollider::inelasticInteractionPossible" << G4endl;
  }

  const G4double coeff = 0.001 * 1.2;
  const G4double one_third = 1.0 / 3.0;

  G4bool possible = true;
  G4double at;
  G4double zt;
  G4double ab;
  G4double zb;

  if(G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {
    at = nuclei_target->getA();
    zt = nuclei_target->getZ(); 
    if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {
      ab = nuclei_bullet->getA();
      zb = nuclei_bullet->getZ();     
    }
    else {

      G4InuclElementaryParticle* particle =
	dynamic_cast<G4InuclElementaryParticle*>(bullet);

      ab = 1;
      zb = particle->getCharge();
    }; 
  }
  else {
    if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) {
      ab = nuclei_bullet->getA();
      zb = nuclei_bullet->getZ();     

      G4InuclElementaryParticle* particle =
	dynamic_cast<G4InuclElementaryParticle*>(target);

      at = 1;
      zt = particle->getCharge();    
    }
    else {

      return possible;
    };  
  }; 

  G4double VCOL = coeff * zt * zb / (pow(at, one_third) + pow(ab, one_third));

  possible = VCOL < ekin;

  return possible;

}
	
G4InteractionCase G4InuclCollider::bulletTargetSetter(G4InuclParticle* bullet,
						      G4InuclParticle* target) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclCollider::bulletTargetSetter" << G4endl;
  }

  G4InteractionCase interCase;

  if(G4InuclNuclei* nuclei_target = dynamic_cast<G4InuclNuclei*>(target)) {     
    if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) { // A + A         
      interCase.setInterCase(2);
      if(nuclei_target->getA() >= nuclei_bullet->getA()) {
	interCase.setBulletTarget(bullet, target);
      }
      else {
	interCase.setBulletTarget(target, bullet);
      }; 
    }
    else {
      interCase.setInterCase(1);
      interCase.setBulletTarget(bullet, target);
    }; 
  }
  else {
    if(G4InuclNuclei* nuclei_bullet = dynamic_cast<G4InuclNuclei*>(bullet)) { 
      if(G4InuclElementaryParticle* part = 
	 dynamic_cast<G4InuclElementaryParticle*>(target)) {
	interCase.setInterCase(1);
	interCase.setBulletTarget(target, bullet);
      };
    };
  };

  return interCase;
}       

G4bool G4InuclCollider::explosion(G4InuclNuclei* target) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4InuclCollider::explosion" << G4endl;
  }

  const G4double a_cut = 20.0;
  const G4double be_cut = 3.0;

  G4double a = target->getA();
  G4double z = target->getZ();
  G4double eexs = target->getExitationEnergy();
  G4bool explo = true;

  if(a > a_cut) {
    explo = false;
  }
  else {
    if(eexs < be_cut * bindingEnergy(a, z)) explo = false;
  };   

  return explo;
}
 
