//#define CHECK_MOMC

#include <iomanip.h>

#include "globals.hh"
//#include "Randomize.hh"

#include "G4Collider.hh"
#include "G4InuclCollider.hh"
#include "G4IntraNucleiCascader.hh"
#include "G4NonEquilibriumEvaporator.hh"
#include "G4EquilibriumEvaporator.hh"
#include "G4Fissioner.hh"
#include "G4BigBanger.hh"
#include "G4ElementaryParticleCollider.hh"
#include "G4InuclParticle.hh"
#include "G4InuclElementaryParticle.hh"
#include "G4InuclNuclei.hh"
#include "G4CollisionOutput.hh"
#include "G4Analyser.hh"
#include "G4WatcherGun.hh"

#include "vector"


G4int testINC();
G4int testINCEvap();
G4int testINCAll();
G4int test();

int main() {

  testINC();     // Only INC model
  testINCEvap(); // INC and evaporation
  //  testINCAll();  // INC, pre-eq, evap, fission 
  test();        // misc. testing

  return 0;       
};

G4int testINCEvap() {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> testINCEvap " << G4endl;
  }
  
  return 0;
};

G4int testINC() {
  G4int verboseLevel = 2; // eguals 1 for data file production, 2 for testing, 3 all

  if (verboseLevel > 1) {
    G4cout << " >>> testINC" << G4endl;
  }

  G4CollisionOutput TRFoutput;

  typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;

  // 0.8 GeV proton with momentum along Z axis
  G4InuclParticle* bullet = new G4InuclElementaryParticle(0.1,1); // momentumBullet, bulletType 

  // Set target
  G4std::vector<G4double> targetMomentum(4, 0.0);

  G4InuclNuclei * target = new G4InuclNuclei(0.,197.,79.); //Au197  momentum = 0
  target->setEnergy();

  // Resigister collider
  G4ElementaryParticleCollider* collider = new G4ElementaryParticleCollider;
  G4IntraNucleiCascader*        cascader = new G4IntraNucleiCascader;
 
  cascader->setElementaryParticleCollider(collider);
  cascader->setInteractionCase(1); // Interaction type is particle with nuclei.

  G4int nCollisions = 10;

    if (verboseLevel > 1) {
  G4cout << setw(6)<< "#ev" << setw(6)  << "part" << setw(11) << "Ekin [GeV]" << setw(11) << "momx" << setw(11) << "momy" << setw(11) << "momz" << G4endl;

    }

  for (G4int i = 1; i <= nCollisions; i++){
    // Make INC
    if (verboseLevel > 2) {
      G4cout << "collision " << i << G4endl; 
    }

    G4CollisionOutput output =  cascader->collide(bullet, target); 

    if (verboseLevel > 2) {
      G4cout << " After Cascade " << G4endl;
      output.printCollisionOutput();
    }
	  
    // Convert Bertini data to Geant4 format
    G4std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();

    if(!particles.empty()) 
      { 
	particleIterator ipart;

	for(ipart = particles.begin(); ipart != particles.end(); ipart++) 
	  {
	    G4std::vector<G4double> mom = ipart->getMomentum();
	    G4double ekin = ipart->getKineticEnergy() * GeV;
	    G4int type = ipart->type();

	    if (verboseLevel > 0) {

    if (verboseLevel > 0) {
  cout.precision(4);

  G4cout << setw(6) << i << setw(6)  << type << setw(11) << ekin << setw(11) << mom[1] * GeV << setw(11) << mom[2] * GeV << setw(11) << mom[3] * GeV << G4endl;
	    }
	  }
      }
    
      }
  }
  return 0;
};


G4int testINCAll()  {

  G4int verboseLevel = 0;

  if (verboseLevel > 1) {
    G4cout << " >>> testINCAll() " << G4endl;
  }
  // Set the default random engine to RanecuEngine
  //  RanecuEngine defaultEngine;
  //HepRandom::setTheEngine(&defaultEngine);
  //  HepRandom::setTheSeed(345354);
  //G4double x = RandFlat::shoot();
  //  G4double koe = G4UniformRand();
  
  // General test program for hetc and inucl

  const G4int to_report = 1;
  G4int nrain = 100; // number of interactions to be generated
  G4double eMin  = 0.5; // minimun energy for bullet
  G4double eMax  = 5.5;   // maximum energy for bullet
  G4int    eBins = 1;   // bullet energy bins
  G4double eStep = (eMax-eMin)/eBins;
  for(G4int e = 0; e < eBins; e++) { // Scan with different energy
 
    // Auxiliarly stuff for ugly analysis
    G4Analyser* analyser = new G4Analyser();
    G4WatcherGun* gun = new G4WatcherGun;
    gun->setWatchers();
    analyser->setWatchers(gun->getWatchers());
    analyser->setInelCsec(1760.0, true);

    // Colliders initialisation
    G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*        cascade = new G4IntraNucleiCascader; // the actual cascade
    G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
    G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
    G4Fissioner*                     fiss = new G4Fissioner;
    G4BigBanger*                     bigb = new G4BigBanger;
    G4InuclCollider*             collider = new G4InuclCollider(colep, cascade, noneq, eqil, fiss, bigb);

    // Bullet / Target initialisation
    // Bullet could be nucleon or pion or nuclei
    // proton momentum in Z-direction [GeV]
    G4double bulletEnergy = eMin + eStep * e; 

    if (verboseLevel > 2) {
      G4cout << "Bullet E =" << bulletEnergy << " GeV" << endl;
    };

    G4InuclParticle* bull = new G4InuclElementaryParticle(bulletEnergy, 1); 

    /*
      Neutron  with the momentum defined by the vector momb(4)
      vector<G4double> momb(4);
      momb[1] = 0.0; 
      momb[2] = 0.0; 
      momb[3] = 2.;
      G4InuclParticle* bull = new G4InuclElementaryParticle(mom, 2);
    */

    /*
      He4 nuclei with momentum momb
      G4InuclParticle* bull = new G4InuclNuclei(momb, 4.0, 2.0);
      ((G4InuclNuclei*)bull)->setEnergy(); // for nuclei mom[4] has to be set specially !!! :::
    */
 
    // Target could be nucleon or nuclei (more precise, in case when
    // particle-particle interaction will be generated, at least one particle
    // has to be nucleon)

    // Au197 target at rest
    // G4InuclParticle* targ = new G4InuclNuclei(0.0, 197.0, 79.0);

    G4InuclParticle* targ = new G4InuclNuclei(0.0, 208.0, 82.0); // Pb

    /*
      // Neutron with momentum momta
      vector<G4double> momta(4);
      momta[1] = -0.3; 
      momta[2] = 0.2; 
      momta[3] = -0.2; 
      G4InuclParticle* targ = new G4InuclElementaryParticle(momta, 2);
    */

    /*
      // C12 nuclei with momentum momta
      G4InuclParticle* targ = new G4InuclNuclei(momta, 12.0, 6.0);
      ((G4InuclNuclei*)targ)->setEnergy();
    */

#ifdef CHECK_MOMC
    vector<G4double> total_mom_in = bull->getMomentum();
    vector<G4double> momt = targ->getMomentum();
    for(G4int i = 0; i < 4; i++) total_mom_in[i] += momt[i];
    vector<G4double> total_mom_out;
    bull->printParticle();
    targ->printParticle();
    if (verboseLevel > 1) {
      G4cout << " tot in mom: px " << total_mom_in[1] << " py " << total_mom_in[2] << " pz " << total_mom_in[3] << " e " << total_mom_in[0] << G4endl;
    }

#endif
    for(G4int i = 0; i < nrain; i++) {
      if((i + 1) % to_report == 0) 

	if (verboseLevel > 2) {
	  G4cout << " Event " << i+1 <<":" << G4endl;
      	}

      G4CollisionOutput cascadeParticles = collider->collide(bull, targ); // standard method
      //  auxiliarly method to get more information about the different stages 
      //  of interaction  
      //  G4CollisionOutput cascadeParticles = collider->collideAndTest(bull, targ, analyser);
      //  G4cout << " ++++++++++++++++++++++++++++++++++++++++++++++++++ " << G4endl;
      //  G4cout << " Global output " << G4endl;
      //  cascadeParticles.printCollisionOutput();
#ifdef CHECK_MOMC
      total_mom_out = cascadeParticles.getTotalOutputMomentum();
      G4cout << " 4 - momentum conservation check " << G4endl
	     << " dE " << total_mom_out[0] - total_mom_in[0] << 
	" dPx " << total_mom_out[1] - total_mom_in[1] <<
	" dPy " << total_mom_out[2] - total_mom_in[2] <<
	" dPz " << total_mom_out[3] - total_mom_in[3] << G4endl;
#endif
      analyser->analyse(cascadeParticles);
    };
    //  analyser->printResults();
    //  analyser->printResultsSimple();
    analyser->printResultsNtuple();
  }

  return 0;
};


G4int test() {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> test() " << G4endl;
  }

  if (verboseLevel > 1) {
    G4cout << " MeV: " << MeV << " GeV: " << GeV << G4endl;
  }

  return 0;
};

