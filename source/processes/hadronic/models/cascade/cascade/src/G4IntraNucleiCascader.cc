//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
#define RUN

#include "G4IntraNucleiCascader.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4LorentzConvertor.hh"
#include "G4ParticleLargerEkin.hh"
#include "G4NucleiModel.hh"
#include "G4CascadParticle.hh"
#include "g4std/algorithm"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

G4IntraNucleiCascader::G4IntraNucleiCascader()
  : verboseLevel(1) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4IntraNucleiCascader::G4IntraNucleiCascader" << G4endl;
  }

}

G4CollisionOutput G4IntraNucleiCascader::collide(G4InuclParticle* bullet,
						 G4InuclParticle* target) {

  if (verboseLevel > 3) {
    G4cout << " >>> G4IntraNucleiCascader::collide" << G4endl;
  }

  const G4int itry_max = 1000;
  const G4int reflection_cut = 500;
  //  const G4double eexs_cut = 0.0001;

  if (verboseLevel > 3) {
    bullet->printParticle();
    target->printParticle();
  }

  G4CollisionOutput output;

#ifdef RUN
  G4InuclNuclei* tnuclei = dynamic_cast<G4InuclNuclei*>(target);
  G4InuclNuclei* bnuclei = dynamic_cast<G4InuclNuclei*>(bullet);
  G4InuclElementaryParticle* bparticle = dynamic_cast<G4InuclElementaryParticle*>(bullet);
  G4NucleiModel model(tnuclei);
  G4std::vector<G4double> momentum_in = bullet->getMomentum();

  momentum_in[0] += tnuclei->getMass();

  G4double ekin_in; 

  if (verboseLevel > 3) {
    model.printModel();
    G4cout << " intitial momentum  E " << momentum_in[0] << " Px " << momentum_in[1] 
	   << " Py " << momentum_in[2] << " Pz " << momentum_in[3] << G4endl;
  }

  G4int itry = 0;

  while (itry < itry_max) {
    itry++;
    model.reset();

    G4std::vector<G4CascadParticle> cascad_particles;
    G4ExitonConfiguration theExitonConfiguration;
    G4std::vector<G4InuclElementaryParticle> output_particles;
    G4double afin = tnuclei->getA();
    G4double zfin = tnuclei->getZ();
   
    if (inter_case == 1) { // particle with nuclei
      ekin_in = bparticle->getKineticEnergy();
      zfin += bparticle->getCharge();
      if (bparticle->nucleon()) afin += 1.0;
      cascad_particles.push_back(model.initializeCascad(bparticle));

    } else { // nuclei with nuclei

      ekin_in = bnuclei->getKineticEnergy();

      G4double ab = bnuclei->getA();
      G4double zb = bnuclei->getZ();

      afin += ab;
      zfin += zb;

      G4std::pair<G4std::vector<G4CascadParticle>, G4std::vector<G4InuclElementaryParticle> > 
	all_particles = model.initializeCascad(bnuclei, tnuclei);

      cascad_particles = all_particles.first;

      for (G4int ip = 0; ip < G4int(all_particles.second.size()); ip++) 
	output_particles.push_back(all_particles.second[ip]);

      if (cascad_particles.size() == 0) { // compound nuclei
	G4int ia = int(ab + 0.5);
	G4int iz = int(zb + 0.5);
	G4int i;

	for (i = 0; i < ia; i++) {
	  G4int knd = i < iz ? 1 : 2;
	  theExitonConfiguration.incrementQP(knd);
	};

	G4int ihn = int(2.0 * (ab - zb) * inuclRndm() + 0.5);
	G4int ihz = int(2.0 * zb * inuclRndm() + 0.5);

	for (i = 0; i < ihn; i++) theExitonConfiguration.incrementHoles(2);

	for (i = 0; i < ihz; i++) theExitonConfiguration.incrementHoles(1);
      };
    }; 

    G4std::vector<G4CascadParticle> new_cascad_particles;
    G4int iloop = 0;

    while (!cascad_particles.empty() && !model.empty()) {
      iloop++;

      if (verboseLevel > 3) {
	G4cout << " Number of cparticles " << cascad_particles.size() << G4endl;
	cascad_particles.back().print();
      }

      new_cascad_particles = model.generateParticleFate(cascad_particles.back(),
							theElementaryParticleCollider);
      if (verboseLevel > 2) {
	G4cout << " ew particles " << new_cascad_particles.size() << G4endl;
      }

      // handle the result of a new step

      if (new_cascad_particles.size() == 1) { // last particle goes without interaction
	cascad_particles.pop_back();

	if (model.stillInside(new_cascad_particles[0])) { // particle survives 

	  if (verboseLevel > 3) {
	    G4cout << " still inside " << G4endl;
	  }

	  if (new_cascad_particles[0].getNumberOfReflections() < reflection_cut &&
	     model.worthToPropagate(new_cascad_particles[0])) { // it's ok

	    if (verboseLevel > 3) {
	      G4cout << " survives " << G4endl;
	    }

	    cascad_particles.push_back(new_cascad_particles[0]);

	  } else { // it becomes an exiton 

	    if (verboseLevel > 3) {
	      G4cout << " becomes an exiton " << G4endl;
	    }

	    theExitonConfiguration.incrementQP(new_cascad_particles[0].getParticle().type());
	  };  

	} else { // goes out

	  if (verboseLevel > 3) {
	    G4cout << " Goes out " << G4endl;
	    new_cascad_particles[0].print();
	  }

	  output_particles.push_back(new_cascad_particles[0].getParticle());
	}; 

      } else { // interaction 

	cascad_particles.pop_back();

	for (G4int i = 0; i < G4int(new_cascad_particles.size()); i++) 
	  cascad_particles.push_back(new_cascad_particles[i]);

	G4std::pair<G4int, G4int> holes = model.getTypesOfNucleonsInvolved();

	theExitonConfiguration.incrementHoles(holes.first);

	if (holes.second > 0) theExitonConfiguration.incrementHoles(holes.second);

      };
    };

    // cascad is finished -> check, whether it's o'k

    if (verboseLevel > 3) {
      G4cout << " Cascade finished  " << G4endl
	     << " output_particles  " << output_particles.size() <<  G4endl;
    }

    G4std::vector<G4double> momentum_out(4, 0.0);
    particleIterator ipart;

    for (ipart = output_particles.begin(); ipart != output_particles.end(); ipart++) {

      G4std::vector<G4double> mom = ipart->getMomentum();

      for (G4int j = 0; j < 4; j++) momentum_out[j] += mom[j];

      zfin -= ipart->getCharge();

      if (ipart->nucleon()) afin -= 1.0;

    };

    if (verboseLevel > 3) {
      G4cout << "  afin " << afin << " zfin " << zfin <<  G4endl;
    }

    if (afin > 1.0) {

      G4InuclNuclei outgoing_nuclei(afin, zfin);
      G4double mass = outgoing_nuclei.getMass();
      momentum_out[0] += mass;        

      for (int j = 0; j < 4; j++) momentum_out[j] = momentum_in[j] - momentum_out[j];

      if (verboseLevel > 3) {
	G4cout << "  Eex + Ekin " << momentum_out[0]  <<  G4endl;
      }

      if (momentum_out[0] > 0.0) { // Eex + Ekin > 0.0
	G4double pnuc = momentum_out[1] * momentum_out[1] + 
	  momentum_out[2] * momentum_out[2] +
	  momentum_out[3] * momentum_out[3]; 
	G4double ekin = sqrt(mass * mass + pnuc) - mass;
	G4double Eex = 1000.0 * (momentum_out[0] - ekin);

	if (verboseLevel > 3) {
	  G4cout << "  Eex  " << Eex  <<  G4endl;
	}

	if (goodCase(afin, zfin, Eex, ekin_in)) { // ok, exitation energy > cut
	  G4std::sort(output_particles.begin(), output_particles.end(), G4ParticleLargerEkin());
	  output.addOutgoingParticles(output_particles);
	  outgoing_nuclei.setMomentum(momentum_out);
	  outgoing_nuclei.setEnergy();
	  outgoing_nuclei.setExitationEnergy(Eex);
	  outgoing_nuclei.setExitonConfiguration(theExitonConfiguration);	                           	  
          output.addTargetFragment(outgoing_nuclei);

	  return output;
	};
      };

    } else { // special case, when one has no nuclei after the cascad

      if (afin == 1.0) { // recoiling nucleon

	for (int j = 0; j < 4; j++) momentum_out[j] = momentum_in[j] - momentum_out[j];

	G4InuclElementaryParticle  last_particle;

	if (zfin == 1.0) { // recoiling proton
	  last_particle.setType(1);

	} else { // neutron

	  last_particle.setType(2);
	}; 

	last_particle.setMomentum(momentum_out);
	output_particles.push_back(last_particle);
      }; 

      G4std::sort(output_particles.begin(), output_particles.end(), G4ParticleLargerEkin());
      output.addOutgoingParticles(output_particles);

      return output;
    }; 
  };

#else

  // special branch to avoid the cascad generation but to get the input for evaporation etc

  G4std::vector<G4double> momentum_out(4, 0.0);
  G4InuclNuclei outgoing_nuclei(169, 69);

  outgoing_nuclei.setMomentum(momentum_out);
  outgoing_nuclei.setEnergy();
  outgoing_nuclei.setExitationEnergy(150.0);

  G4ExitonConfiguration theExitonConfiguration(3.0, 3.0, 5.0, 6.0);

  outgoing_nuclei.setExitonConfiguration(theExitonConfiguration);	                           
  output.addTargetFragment(outgoing_nuclei);

  return output;

  /*
    G4InuclElementaryParticle* bparticle = dynamic_cast<G4InuclElementaryParticle*>
    (bullet);
    G4InuclNuclei* tnuclei = dynamic_cast<G4InuclNuclei*>(target);
    output.addOutgoingParticle(*bparticle);
    output.addTargetFragment(*tnuclei);
    return output;
  */

#endif

  if (verboseLevel > 3) {
    G4cout << " IntraNucleiCascader-> no inelastic interaction after " << itry_max << " attempts "
	   << G4endl;
  }

  output.trivialise(bullet, target);

  return output;
}

G4bool G4IntraNucleiCascader::goodCase(G4double a, 
				       G4double z, 
				       G4double eexs, 
				       G4double ein) const {

  if (verboseLevel > 3) {
    G4cout << " >>> G4IntraNucleiCascader::goodCase" << G4endl;
  }

  const G4double eexs_cut = 0.0001;
  const G4double reason_cut = 7.0;
  const G4double ediv_cut = 5.0;

  G4bool good = false;

  if (eexs > eexs_cut) {
    G4double eexs_max0z = 1000.0 * ein / ediv_cut;
    G4double dm = bindingEnergy(a, z);
    G4double eexs_max = eexs_max0z > reason_cut*dm ? eexs_max0z : reason_cut * dm;

    if(eexs < eexs_max) good = true;

    if (verboseLevel > 3) {
      G4cout << " eexs " << eexs << " max " << eexs_max << " dm " << dm << G4endl;
    }

  };

  return good; 
}
