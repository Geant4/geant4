
#include "globals.hh"
#include "G4IonTable.hh"

#include "G4InuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4InuclParticle.hh"
#include "G4CollisionOutput.hh"
#include "G4Nucleus.hh"
#include "G4NucleiModel.hh"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef G4std::vector<G4InuclNuclei>::iterator nucleiIterator;

int main(int argc, char **argv ) {

  G4int verboseLevel = 3;

  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, photon = 10 };

  G4int bulletType = 0;

  G4std::vector<G4double>  momentumBullet(4, 0.0);
  momentumBullet[0] = 1.37126;
  momentumBullet[3] = 1.5;

  G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, 1); 

  G4InuclNuclei*   target  = NULL;
  G4InuclParticle* targetH = NULL;
  G4NucleiModel*   model = NULL;

  G4double theNucleusA = 1;
  G4std::vector<G4double> targetMomentum(4, 0.0);

  G4CollisionOutput output;

  G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
  G4IntraNucleiCascader*            inc = new G4IntraNucleiCascader; // the actual cascade
  inc->setInteractionCase(1); // Interaction type is particle with nuclei.

  G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
  G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
  G4Fissioner*                     fiss = new G4Fissioner;
  G4BigBanger*                     bigb = new G4BigBanger;

  G4InuclCollider*             collider = new G4InuclCollider(colep, inc, noneq, eqil, fiss, bigb);

  for (G4int i = 1; i< 100 ; i++) {
    if ( theNucleusA < 1.5 ) 
      {
	model = new G4NucleiModel(new G4InuclNuclei(targetMomentum, 1, 1));
	targetH = new G4InuclElementaryParticle((model->generateNucleon(1, 1)).getMomentum(), 1); 
   
	//		do
	  {
	    cout << "+";
	    output = collider->collide(bullet, targetH); 
	  } 
	  //		while(output.getOutgoingParticles().size()<2.5);
      } 
    else 
      {
	output = collider->collide(bullet, target ); 
      }

    if (verboseLevel > 1) 
      {
	G4cout << " Cascade output: " << G4endl;
	output.printCollisionOutput();
      }
  
    // Convert cascade data to use hadronics interface

    G4std::vector<G4InuclNuclei>             nucleiFragments = output.getNucleiFragments();
    G4std::vector<G4InuclElementaryParticle> particles =       output.getOutgoingParticles();

    G4int numSecondaries = nucleiFragments.size()+particles.size();
    cout << "num secondaries: " << numSecondaries << G4endl;
    if(!particles.empty()) { 
      particleIterator ipart;
      G4int outgoingParticle;

      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	outgoingParticle = ipart->type();
	G4std::vector<G4double> mom = ipart->getMomentum();
	G4double ekin = ipart->getKineticEnergy() * GeV;
	G4ThreeVector aMom(mom[1], mom[2], mom[3]);
	aMom = aMom.unit();

      
      }

      nucleiIterator ifrag;

      for(ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) 
	{
	  G4double eKin = ifrag->getKineticEnergy() * GeV;
	  G4std::vector<G4double> mom = ifrag->getMomentum();
	  G4ThreeVector aMom(mom[1], mom[2], mom[3]);
	  aMom = aMom.unit();

	  // hpw @@@ ==> Should be zero: G4double fragmentExitation = ifrag->getExitationEnergyInGeV();

	  if (verboseLevel > 2) {
	    G4cout << " Nuclei fragment: " << G4endl;
	    ifrag->printParticle();
	  }
	  G4int A = G4int(ifrag->getA());
	  G4int Z = G4int(ifrag->getZ());

	}
    }
  }
}




