//#define CHECK_MOMC

#include <iomanip.h>

#include "globals.hh"
#include "Randomize.hh"

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
#include "G4ThreeVector.hh"

typedef G4std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef G4std::vector<G4InuclNuclei>::iterator nucleiIterator;

enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

G4int nCollisions = 10;      // collisions to be generated
G4int  bulletType = proton;    // bullet particle
G4double     momZ = 160;      // momentum in z-direction
G4double        A = 27.0;      // target atomic weight Al
G4double        Z = 13.0;      // target atomic number

G4CollisionOutput output;

G4int testINC(G4int, G4int, G4double, G4double, G4double);
G4int testINCEvap();
G4int testINCAll(G4int, G4int, G4double, G4double, G4double);
G4int printData(G4int event);
G4int test();

int main(int argc, char **argv ) {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> cascade::main " << G4endl;
  }

  // Get argumets from command line
  nCollisions =           (argc > 1) ? atoi(argv[1]) : nCollisions;
  bulletType  =           (argc > 2) ? atoi(argv[2]) : proton;
  momZ        = G4double(((argc > 3) ? atoi(argv[3]) : momZ)) / GeV;
  A           = G4double(((argc > 4) ? atoi(argv[4]) : A));
  Z           = G4double(((argc > 5) ? atoi(argv[5]) : Z));

  if (verboseLevel > 1) {
    G4cout << " nCollisions " << nCollisions << G4endl;
    G4cout << "  bulletType " << bulletType  << G4endl;
    G4cout << "        momZ " << momZ        << G4endl;
    G4cout << "           A " << A           << G4endl;
    G4cout << "           Z " << Z           << G4endl;
  }

  //testINC(nCollisions, bulletType, momZ, A, Z);     // Only INC model  
  testINCAll(nCollisions, bulletType, momZ, A, Z);  // INC, pre-eq, evap, fission
 
  //testINCEvap(); // INC and evaporation 
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

G4int testINC(G4int nCollisions, G4int bulletType, G4double momZ, G4double A, G4double Z) {
  G4int verboseLevel = 1; // eguals 1 for data file production, 2 for testing, 3 all

  if (verboseLevel > 1) {
    G4cout << " >>> testINC" << G4endl;
  }

  G4CollisionOutput TRFoutput;

  G4InuclParticle* bullet = new G4InuclElementaryParticle(momZ, bulletType); // momentumBullet, bulletType 

  // Set target
  G4std::vector<G4double> targetMomentum(4, 0.0);

    G4InuclNuclei*   target  = NULL;
    G4InuclParticle* targIsH = NULL;

    if ( A == 1 ) {
      targIsH = new G4InuclElementaryParticle(0.0, 1);     
    } else { 
      target = new G4InuclNuclei(0.0, A, Z);
      target->setEnergy();
    }

  // Resigister collider
  G4ElementaryParticleCollider* collider = new G4ElementaryParticleCollider;
  G4IntraNucleiCascader*        cascader = new G4IntraNucleiCascader;
 
  cascader->setElementaryParticleCollider(collider);
  cascader->setInteractionCase(1); // Interaction type is particle with nuclei.

  if (verboseLevel > 1) {
    G4cout << 
      setw(6)  << "#ev"        << 
      setw(6)  << "part"       << 
      setw(11) << "Ekin [GeV]" << 
      setw(11) << "momx"       << 
      setw(11) << "momy"       << 
      setw(11) << "momz"       << G4endl;
  }

  for (G4int i = 1; i <= nCollisions; i++) {

    if (verboseLevel > 0) {
      G4cout << "collision " << i << G4endl; 
    }

    if ( A == 1 ) {
      output = cascader->collide(bullet, targIsH); // standard method
    } else {
      output = cascader->collide(bullet, target); // Make INC    
    }
    
    printData(i);
  }

  return 0;
};

G4int testINCAll(G4int nCollisions, G4int bulletType, G4double momZ, G4double A, G4double Z) {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> testINCAll() " << G4endl;
  }

  G4double eMin  = 0.5;   // minimun energy for bullet
  G4double eMax  = 5.5;   // maximum energy for bullet
  G4int    eBins = 1;     // bullet energy bins
  G4double eStep = (eMax - eMin) / eBins;

  for(G4int e = 0; e < eBins; e++) { // Scan with different energy
 
    // Colliders initialisation
    G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*        cascade = new G4IntraNucleiCascader; // the actual cascade
    cascade->setInteractionCase(1); // Interaction type is particle with nuclei.

    G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
    G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
    G4Fissioner*                     fiss = new G4Fissioner;
    G4BigBanger*                     bigb = new G4BigBanger;
    G4InuclCollider*             collider = new G4InuclCollider(colep, cascade, noneq, eqil, fiss, bigb);

    // Bullet / Target initialisation
    G4double bulletEnergy = eMin + eStep * e; 

    if (verboseLevel > 2) {
      G4cout << "Bullet E =" << bulletEnergy << " GeV" << endl;
    };

    // Set target
    G4std::vector<G4double> targetMomentum(4, 0.0);

    G4InuclParticle* bull = new G4InuclElementaryParticle(momZ, proton);  

    G4InuclNuclei* targ = NULL;
    G4InuclParticle* targIsH = NULL;

    if ( A == 1 ) {
      targIsH = new G4InuclElementaryParticle(0.0, 2.0);     
    } else { 
      targ = new G4InuclNuclei(0.0, A, Z);
      targ->setEnergy();      
    }

    if (verboseLevel > 2) {
      G4cout << " Event " << e+1 <<":" << G4endl;
    }

    for (G4int i = 1; i <= nCollisions; i++) {

      if (verboseLevel > 3) {
	G4cout << "collision " << i << G4endl; 
      }

    if ( A == 1 ) {
      output = collider->collide(bull, targIsH); // standard method
    } else {
      output = collider->collide(bull, targ); // standard method
    }
      printData(i);
    }

    delete bull;
    delete targ;

    delete colep;
    delete cascade; 
    delete noneq;
    delete eqil;
    delete fiss;
    delete bigb;
    delete collider;
  }

  return 0;
};

G4int printData(G4int i) {

  G4int verboseLevel = 1;

  if (verboseLevel > 2) {
    G4cout << " After Cascade " << G4endl;
    output.printCollisionOutput();
  }
	  
  // Convert Bertini data to Geant4 format

  G4std::vector<G4InuclNuclei> nucleiFragments = output.getNucleiFragments();

  if(!nucleiFragments.empty()) { 
    nucleiIterator ifrag;
    
    
    for(ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) {
      G4std::vector<G4double> m = ifrag->getMomentum();
      G4ThreeVector mom(m[1], m[2], m[3]);    
      G4double ekin = ifrag->getKineticEnergy() * GeV;

      G4int type = 0; // :::

      if (verboseLevel > 2) {
	G4cout << " Fragment mass: " << ifrag->getMass()  << G4endl;
        G4cout << " Momentum magnitude: " << mom.mag() << G4endl;
      }

      G4double fEx = ifrag->getExitationEnergyInGeV();
      G4int fA = G4int(ifrag->getA());
      G4int fZ = G4int(ifrag->getZ());

      if (verboseLevel > 2) {

	G4cout << " Nuclei fragment: " << G4endl;
	G4cout << " exitation energy " << fEx << " " << fA << " " << fZ << G4endl;
 
        ifrag->printParticle();
      }
	
      if (verboseLevel > 0) {
	cout.precision(4);

	G4cout << 
	  setw(6)  << i            << 
	  setw(6)  << type         << 
	  setw(11) << ekin         << 
	  setw(11) << mom[0] * GeV << 
	  setw(11) << mom[1] * GeV << 
	  setw(11) << mom[2] * GeV << 
	  setw(13) << fA           << 
	  setw(13) << fZ           << 
	  setw(13) << fEx          << G4endl;
      }
    }
  }

  G4std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();
  if(!particles.empty()) { 
    particleIterator ipart;

    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
      G4std::vector<G4double> mom = ipart->getMomentum();
      G4double ekin = ipart->getKineticEnergy() * GeV;
      G4int type = ipart->type();

      if (verboseLevel > 0) {
	cout.precision(4);

	G4cout << 
	  setw(6)  << i            << 
	  setw(6)  << type         << 
	  setw(11) << ekin         << 
	  setw(11) << mom[1] * GeV << 
	  setw(11) << mom[2] * GeV << 
	  setw(11) << mom[3] * GeV << 
	  setw(13) << 0            << 
	  setw(13) << 0            << 
	  setw(13) << 0.0          << G4endl;
      }
    }
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






