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
//#define CHECK_MOMC

#include <vector>
#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

#include "globals.hh"
#include "Randomize.hh"

#include "G4Collider.hh"
#include "G4InuclCollider.hh"
#include "G4PreCompoundInuclCollider.hh"
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
#include "G4ThreeVector.hh"
#include "G4NucleiModel.hh"
#include "G4LorentzConvertor.hh"

#include "G4InuclCollider.hh"
#include "G4PreCompoundInuclCollider.hh"

G4bool coulombOK;
G4double kE=1000; // scale energy ok

G4double eInit = 0.0;
G4double eTot = 0.0;
G4double sumBaryon = 0.0;
G4double sumEnergy = 0.0;


typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;

enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

G4int runId       = 0; 
G4int nCollisions = 10;      // collisions to be generated
G4int  bulletType = proton;    // bullet particle
G4double     momZ = 160;      // momentum in z-direction [MeV/c]
G4double        A = 27.0;      // target atomic weight Al
G4double        Z = 13.0;      // target atomic number

G4int n1 = 0;
G4int n2 = 0;
G4int n3 = 0;
G4int n5 = 0;
G4int n7 = 0;
G4int n10 = 0;

G4CollisionOutput output;

G4int tCoulomb(G4int, G4int, G4int, G4double, G4double, G4double);
//G4int testINC(G4int, G4int, G4double, G4double, G4double);
G4int testINCEvap();
G4int testINCAll(G4int, G4int, G4double, G4double, G4double);
G4int printData(G4int runId, G4int event, G4int cok);
G4int printCross(G4int i);

G4int test();

G4int verboseLevel = 2;
G4InuclElementaryParticle* bull;

int main(int argc, char **argv ) {

  if (verboseLevel > 3) {
    G4cout << " >>> cascade::main " << G4endl;
  }

  // Get argumets from command line
  runId       =           (argc > 1) ? atoi(argv[1]) : runId;
  nCollisions =           (argc > 2) ? atoi(argv[2]) : nCollisions;
  bulletType  =           (argc > 3) ? atoi(argv[3]) : proton;
  momZ        = G4double(((argc > 4) ? atoi(argv[4]) : momZ)) / GeV;  // MeV = 1 an GeV = 1000
  A           = G4double(((argc > 5) ? atoi(argv[5]) : A));
  Z           = G4double(((argc > 6) ? atoi(argv[6]) : Z));

  if (verboseLevel > 3) {
    G4cout << " # collisions " << nCollisions << G4endl;
    G4cout << "  bullet type " << bulletType  << G4endl;
    G4cout << "     momentum " << momZ*kE    << " [MeV]" << G4endl; // kE = 1000; 
    G4cout << "            A " << A           << G4endl;
    G4cout << "            Z " << Z           << G4endl;
  }

  //    G4cout << " MeV: " << MeV << " GeV: " << GeV << G4endl;

  if (verboseLevel > 3) { // deside if heade will be written
    //    G4cout << "# cascade.cc with parameters : runId, nCollisions, bulletType,  momZ,  targetA, targetZ" << G4endl;
    G4cout << runId << " " << nCollisions << " " << bulletType  << " " << momZ*kE  << " " << A << " " << Z << G4endl;
  }

  tCoulomb(runId, nCollisions, bulletType, momZ, A, Z);  // test coulomb
  //testINC(nCollisions, bulletType, momZ, A, Z);     // Only INC model  
  //  testINCAll(nCollisions, bulletType, momZ, A, Z);  // INC, pre-eq, evap, fission
 
  //testINCEvap(); // INC and evaporation 
  //  test();        // misc. testing

  return 0;       
}

G4int testINCEvap() {

  G4int verboseLevel = 1;

  if (verboseLevel > 2) {
    G4cout << " >>> testINCEvap " << G4endl;
  }
  
  return 0;
}

G4int tCoulomb(G4int runId, G4int nCollisions, G4int bulletType, G4double momZ, G4double A, G4double Z) {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> tCoulomb  Start" << G4endl;
  }
  
    n1 = 0;
    n2 = 0;
    n3 = 0;
    n5 = 0;
    n7 = 0;
    n10 = 0;
 
    // Colliders initialisation
    G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*        cascade = new G4IntraNucleiCascader; // the actual cascade
    cascade->setInteractionCase(1); // Interaction type is particle with nuclei.

    G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
    G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
   G4Fissioner*                     fiss = new G4Fissioner;
    G4BigBanger*                     bigb = new G4BigBanger;
    G4InuclCollider*             collider = new G4InuclCollider(colep, cascade, noneq, eqil, fiss, bigb);

    // G4PreCompoundInuclCollider*  collider = new G4PreCompoundInuclCollider(colep, cascade, noneq, bigb);

    std::vector<G4double> targetMomentum(4, 0.0);
    //    std::vector<G4double>  bulletMomentum(4, 0.0);
    //    G4double mass = 0.93827;

    // Old
    //    bulletMomentum[3] = momZ;
    //bulletMomentum[3] = std::sqrt(bulletMomentum[3] * bulletMomentum[3] + 2 * bulletMomentum[3] * mass); // only this is used in tests
    //bulletMomentum[2] = std::sqrt(bulletMomentum[2] * bulletMomentum[2] + 2 * bulletMomentum[2] * mass);
    //bulletMomentum[1] = std::sqrt(bulletMomentum[1] * bulletMomentum[1] + 2 * bulletMomentum[1] * mass); 
  G4double px,py,pz;
  px=0.0 ;
  py=0.0 ;
  pz=momZ/2; //[GeV]

  std::vector<G4double> momentumBullet(4);
  momentumBullet[0] =0.;
  momentumBullet[1] =0;
  momentumBullet[2] =0;
  momentumBullet[3] =std::sqrt(px*px+py*py+pz*pz);

  //  G4cout << "momentumBullet[3]" << momentumBullet[3] << G4endl;
  bull = new G4InuclElementaryParticle(momentumBullet, bulletType); // counts mom[0] = E tot from mom[1]-mom[3]
   


    if (verboseLevel > 2) {
      G4cout << "Bullet:  " << G4endl;  
      bull->printParticle();
    }

    G4InuclNuclei* targ = NULL;
    G4InuclParticle* targIsH = NULL;

    if ( !(G4int(A) == 1) ) {
      targ = new G4InuclNuclei(targetMomentum, A, Z);
      targ->setEnergy();      

      std::vector<G4double>  bmom = bull->getMomentum();
      eInit = std::sqrt(bmom[0] * bmom[0]);
      std::vector<G4double> tmom = targ->getMomentum();
      eInit += std::sqrt(tmom[0] * tmom[0]);

      if (verboseLevel > 2) {
	G4cout << "Target:  " << G4endl;  
	targ->printParticle();
      }
    };

    G4int i;
    for (i = 1; i <= nCollisions; i++) {
      
      if (verboseLevel > 3) {
	G4cout << "collision " << i << G4endl; 
      }

      if ( G4int(A) == 1 ) {
	G4int is = 0;
	targIsH = new G4InuclElementaryParticle(targetMomentum, 1);

	std::vector<G4double>  bmom = bull->getMomentum();
	eInit = std::sqrt(bmom[0] * bmom[0]);
	std::vector<G4double> tmom = targIsH->getMomentum();
	eInit += std::sqrt(tmom[0] * tmom[0]);

	do {
	  if (verboseLevel > 1) {
	    G4cout << "Target:  " << G4endl;  
	    targIsH->printParticle();
	  }

	  output = collider->collide(bull, targIsH);
	  is = output.getOutgoingParticles().size();
	  //	  		cout << "+";
	
	} while( is == 2);
	//		cout << " "  << G4endl;
      } else {


	//	output = collider->collide(bull, targ); 
	//	int maxTries=10;
	int nTries=0;

	
	//    do  // we try to create inelastic interaction 
	  //{
	//	output = collider->collide(bullet, target );
 coulombOK= true;
	output = collider->collide(bull, targ);

	//-------------------------------------
	// test if protn energy too low

   std::vector<G4InuclElementaryParticle> p= output.getOutgoingParticles();
  if(!p.empty()) { 
    for(    particleIterator ipart = p.begin(); ipart != p.end(); ipart++) {
      if (ipart->type() == proton) {
	G4double e = ipart->getKineticEnergy();
      if (e < 0.005){
	//	           G4cout << std::setw(8)  << e    << G4endl; 
	coulombOK= false;
      };

      
      }
    }
  }
  // G4cout << nTries << G4endl;
  //  G4cout << "coulomb: "<<coulombOK << G4endl;
    //-----------------------
	nTries++;

	//   } while (
	       //	      (nTries < maxTries)                                                               &&
	       //(output.getOutgoingParticles().size() + output.getNucleiFragments().size() < 2.5) &&
	       //(output.getOutgoingParticles().size()!=0)                                         &&
	       //(output.getOutgoingParticles().begin()->type()==bullet->type())
	       //);

	       //	       (nTries < maxTries) &&
	       //output.getOutgoingParticles().size() == 1 &&     // we retry when elastic collision happened
               //output.getNucleiFragments().size() == 1 &&            
	       //output.getOutgoingParticles().begin()->type() == bull->type() &&
	       //output.getNucleiFragments().begin()->getA() == targ->getA() && 
	       //output.getNucleiFragments().begin()->getZ() == targ->getZ() &&
	//     coulombOK
	//       );



      }
      //  G4cout << ">>>>coulomb: "<<coulombOK << G4endl;
  printData(runId, i, (int)coulombOK);
      //      printCross(i);
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

    G4double nc = G4double(nCollisions);

    if (verboseLevel > 3) {
      G4cout << 
	std::setw(8)  << momZ*kE    << 
	std::setw(8)  << i    << 
	std::setw(8)  << n1 / nc   << 
	std::setw(13) << n2 / nc   << 
	std::setw(13) << n3 / nc   << 
	std::setw(13) << n5 / nc   << 
	std::setw(13) << n7 / nc   << 
	std::setw(13) << n10 / nc  << G4endl;
    }
  return 0;
}

G4int testINCAll(G4int nCollisions, G4int bulletType, G4double momZ, G4double A, G4double Z) {

  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> testINCAll() " << G4endl;
  }

  //G4double eMin  = 0.5;   // minimun energy for bullet
  //G4double eMax  = 5.5;   // maximum energy for bullet
  G4int    eBins = 1;     // bullet energy bins
  //  G4double eStep = (eMax - eMin) / eBins;

  for (G4int e = 0; e < eBins; e++) { // Scan with different energy

    n1 = 0;
    n2 = 0;
    n3 = 0;
    n5 = 0;
    n7 = 0;
    n10 = 0;
 
    // Colliders initialisation
    G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
    G4IntraNucleiCascader*        cascade = new G4IntraNucleiCascader; // the actual cascade
    cascade->setInteractionCase(1); // Interaction type is particle with nuclei.

    G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
    G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
    G4Fissioner*                     fiss = new G4Fissioner;
    G4BigBanger*                     bigb = new G4BigBanger;
    G4InuclCollider*             collider = new G4InuclCollider(colep, cascade, noneq, eqil, fiss, bigb);

    //    // Bullet / Target initialisation
    //G4double bulletEnergy = eMin + eStep * e; 

    //    if (verboseLevel > 2) {
    //  G4cout << "Bullet E =" << bulletEnergy << " GeV" << G4endl;
    //};

    // Set target
    std::vector<G4double> targetMomentum(4, 0.0);
    std::vector<G4double>  bulletMomentum(4, 0.0);
    G4double mass = 0.93827;
    bulletMomentum[3] = momZ;
    bulletMomentum[3] = std::sqrt(bulletMomentum[3] * bulletMomentum[3] + 2 * bulletMomentum[3] * mass); // only this is used in tests
    bulletMomentum[2] = std::sqrt(bulletMomentum[2] * bulletMomentum[2] + 2 * bulletMomentum[2] * mass);
    bulletMomentum[1] = std::sqrt(bulletMomentum[1] * bulletMomentum[1] + 2 * bulletMomentum[1] * mass); 

    bull = new G4InuclElementaryParticle(bulletMomentum, bulletType); // counts mom[0] = E tot from mom[1]-mom[3]
   
    if (verboseLevel > 2) {
      G4cout << "Bullet:  " << G4endl;  
      bull->printParticle();
    }

    G4InuclNuclei* targ = NULL;
    G4InuclParticle* targIsH = NULL;

    if ( !(G4int(A) == 1) ) {
      targ = new G4InuclNuclei(targetMomentum, A, Z);
      targ->setEnergy();      

      std::vector<G4double>  bmom = bull->getMomentum();
      eInit = std::sqrt(bmom[0] * bmom[0]);
      std::vector<G4double> tmom = targ->getMomentum();
      eInit += std::sqrt(tmom[0] * tmom[0]);

      if (verboseLevel > 2) {
	G4cout << "Target:  " << G4endl;  
	targ->printParticle();
      }

    };

    if (verboseLevel > 2) {
      G4cout << " Event " << e+1 <<":" << G4endl;
    }

    G4int i;
    for (i = 1; i <= nCollisions; i++) {
      
      if (verboseLevel > 4) {
	G4cout << "collision " << i << G4endl; 
      }

      if ( G4int(A) == 1 ) {
	G4int is = 0;
	targIsH = new G4InuclElementaryParticle(targetMomentum, 1);

	std::vector<G4double>  bmom = bull->getMomentum();
	eInit = std::sqrt(bmom[0] * bmom[0]);
	std::vector<G4double> tmom = targIsH->getMomentum();
	eInit += std::sqrt(tmom[0] * tmom[0]);

	do {
	  if (verboseLevel > 1) {
	    G4cout << "Target:  " << G4endl;  
	    targIsH->printParticle();
	  }

	  output = collider->collide(bull, targIsH);
	  is = output.getOutgoingParticles().size();
	  //	  		cout << "+";
	
	} while( is == 2);
	//		cout << " "  << G4endl;
      } else {
	output = collider->collide(bull, targ); 
      }
      //printData(i);
      printCross(i);
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

    G4double nc = G4double(nCollisions);

    if (verboseLevel > -1) {
      G4cout << 
	std::setw(8)  << momZ*kE    << 
	std::setw(8)  << i    << 
	std::setw(8)  << n1 / nc   << 
	std::setw(13) << n2 / nc   << 
	std::setw(13) << n3 / nc   << 
	std::setw(13) << n5 / nc   << 
	std::setw(13) << n7 / nc   << 
	std::setw(13) << n10 / nc  << G4endl;
    }
  }

  return 0;
}

G4int printData(G4int runId, G4int i, G4int cok) {

  G4int verboseLevel = 1;
  if (verboseLevel > 4) {
    G4cout << " After Cascade " << G4endl;
    output.printCollisionOutput();
  }

  sumEnergy = bull->getKineticEnergy(); // In GeV 
  sumBaryon = A;
  if (bulletType == proton || bulletType == neutron) {
    sumBaryon += 1;
  } 

  // Convert Bertini data to Geant4 format
  std::vector<G4InuclNuclei> nucleiFragments = output.getNucleiFragments();

  G4double eKinTot = 0.0;
  G4double ekin = 0.0;
 
  if(!nucleiFragments.empty()) { 
    nucleiIterator ifrag;
    eTot = 0;
    for(ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) {
    
      std::vector<G4double> m = ifrag->getMomentum();

      eTot  += std::sqrt(m[0] * m[0]);

      G4ThreeVector mom(m[1], m[2], m[3]);    
      ekin = ifrag->getKineticEnergy();

      G4int type = 0; 

      if (verboseLevel > 2) {
	G4cout << " Fragment mass: " << ifrag->getMass()  << G4endl;
        G4cout << " Momentum magnitude: " << mom.mag() << G4endl;
      }

      G4double fEx = ifrag->getExitationEnergyInGeV();
      G4int fA = G4int(ifrag->getA());
      G4int fZ = G4int(ifrag->getZ());
      G4int modelId = ifrag->getModel();
      sumBaryon -= fA;
      sumEnergy -= ekin;
      sumEnergy -= fEx;

      if (verboseLevel > 2) {

	G4cout << " Nuclei fragment: " << G4endl;
	G4cout << " exitation energy " << fEx << " " << fA << " " << fZ << G4endl;
 
        ifrag->printParticle();
      }
	
      if (verboseLevel > 0) {
	G4cout.precision(3);

	G4cout << 
	  std::setw(8)  << runId        << 
	  std::setw(8)  << i            << 
	  std::setw(8)  << type         << 
	  std::setw(8)  << modelId      <<
	  std::setw(13) << ekin *kE     << 
	  std::setw(13) << mom[0] *kE   << 
	  std::setw(13) << mom[1] *kE   << 
	  std::setw(13) << mom[2] *kE   << 
	  std::setw(13) << fA           << 
	  std::setw(13) << fZ           << 
	  std::setw(13) << fEx    *kE   <<
	  " "  << cok   << G4endl;
	//	  std::setw(13) << sumBaryon    << 
	//std::setw(13) << sumEnergy    << G4endl;
      }

      eKinTot += ekin;
    }
  }

  std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();
  if(!particles.empty()) { 
    particleIterator ipart;

    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
      std::vector<G4double> mom = ipart->getMomentum();
      eTot   += std::sqrt(mom[0] * mom[0]);

      // std::vector<G4double>  mom(3, 0.0);
      ekin = ipart->getKineticEnergy();
      G4int type = ipart->type();
      G4int modelId = ipart->getModel();
      if (type == proton || type == neutron) {
	sumBaryon -= 1;
      } 

      sumEnergy -= ekin;

      if (verboseLevel > 0) {
	G4cout.precision(4);

	G4cout << 
	  std::setw(8)  << runId        << 
	  std::setw(8)  << i            << 
	  std::setw(8)  << type         <<
	  std::setw(8)  << modelId     <<
	  std::setw(13) << ekin   *kE << 
       	  std::setw(13) << mom[1] *kE << 
       	  std::setw(13) << mom[2] *kE << 
	  std::setw(13) << mom[3] *kE << 
	  std::setw(13) << 0            << 
	  std::setw(13) << 0            << 
	  std::setw(13) << 0.0          << 
	  " "  << cok   << G4endl;
	//	  std::setw(13) << sumBaryon    << 
	// std::setw(13) << sumEnergy    << G4endl;
      }
      eKinTot += ekin;
    }
  }

  if (sumBaryon != 0) {
    G4cout << "ERROR: no baryon number conservation, sum of baryons = " << sumBaryon << G4endl;
  }

  if (verboseLevel > 2) {
    if (sumEnergy > 0.01 ) {
      G4cout << "NOTE: Kinetic energy conservation violated by " << sumEnergy << " GeV" << G4endl;
    }

    G4cout << "ERROR: nergy conservation at level  ~" << (eInit  - eTot) * GeV << " MeV" << G4endl;
  }

  if (sumEnergy < -5.0e-5 ) { // 0.05 MeV
    //    G4cout << "FATAL ERROR: energy created  " << sumEnergy * GeV << " MeV" << G4endl; // this really givesignal ::: 
  }

  if (verboseLevel > 2) {

    G4cout << "Total energy           : " << eTot*kE << "MeV" << G4endl; 
    G4cout << "Total kinetice energy  : " << eKinTot*kE << "MeV" << G4endl; 
    G4cout << "Baryon sum             : " << sumBaryon << G4endl;
  }
  //  eTot = 0.0;
  return 0;
}

G4int printCross(G4int i) {

  G4int verboseLevel = 1;

  if (verboseLevel > 4) {

    G4cout << " After Cascade " << G4endl;
    G4cout << " i " << i << G4endl;
    output.printCollisionOutput();
  }

  // Convert Bertini data to Geant4 format
  std::vector<G4InuclNuclei> nucleiFragments = output.getNucleiFragments();

  if(!nucleiFragments.empty()) { 
    nucleiIterator ifrag;
        
    for(ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) {    
      std::vector<G4double> m = ifrag->getMomentum();
    }
  }

  std::vector<G4InuclElementaryParticle> particles = output.getOutgoingParticles();
  if(!particles.empty()) { 
    particleIterator ipart;
    G4int type;

    for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
      type = ipart->type();

      switch(type) {
      case proton:
	n1 +=1;
	break;
      case neutron:
	n2 +=1;
	break;
      case pionPlus:
	n3 +=1;
	break;
      case pionMinus:
	n5 +=1;
	break;
      case pionZero:
	n7 +=1;
	break;
      case foton:
	n10 +=1;
	break;
      default:
	G4cout << "ERROR: unknown particle" << G4endl;
      }
    }
  }
  return 0;
    
}

G4int test() {
  G4int verboseLevel = 1;

  if (verboseLevel > 2) {
    G4cout << " >>> test() " << G4endl;
  }

  if (verboseLevel > 2) {
    G4cout << " MeV: " << MeV << " GeV: " << GeV << G4endl;
  }

  if (verboseLevel > 2) {
    std::vector<G4double>  m(4, 0.0);

    G4double mZ = 0.585;
    G4double mY = 0.0;
    G4double mX = 0.0;
    G4double mass = 0.93827;

    G4double e = std::sqrt(mZ * mZ + mY * mY + mX * mX + mass * mass);

    m[3] = mZ;

    //      G4double ekin = ipart->getKineticEnergy() * GeV;
    //G4ThreeVector aMom(mom[1], mom[2], mom[3]);
    //aMom = aMom.unit();
  
    G4cout << G4endl << ">>> previous bug in kin e" << G4endl;
    bull = new G4InuclElementaryParticle(m, 1);
    bull->printParticle();

    G4cout << G4endl << ">>> kinetic energy ok in z-dir" << G4endl;
    m[3] = std::sqrt(momZ * momZ + 2 * momZ * mass);
    bull = new G4InuclElementaryParticle(m, 1);
    bull->printParticle();

    G4cout << G4endl << ">>> hole vectos set" << G4endl;

    m[3] = mZ;
    m[2] = 0;
    m[1] = 0;
    m[0] = std::sqrt(m[1] * m[1] + m[2] * m[2] + m[3] * m[3] + mass * mass);

    // fix m so that  ekin with the mass gets ok.
    m[3] = mZ;
    m[2] = mY;
    m[1] = mX;
    m[0] = e;

    G4cout << G4endl << ">>> inciming:" << G4endl;
    bull = new G4InuclElementaryParticle(m, 1); // expects full mom[0]-mom[3] with correct E tot
    bull->printParticle();

    G4double pLength = std::sqrt(m[0] * m[0] - mass * mass); G4cout << " pLength " << pLength  << G4endl;
    G4double mLength = std::sqrt(m[1] * m[1] + m[2] * m[2] + m[3] * m[3]); G4cout << " mLength " << mLength << G4endl;

    G4double scale = 1;
    m[3] = m[3] * scale;
    m[2] = m[2] * scale;
    m[1] = m[1] * scale;

    m[3] = std::sqrt(m[3] * m[3] + 2 * m[3] * mass);
    m[2] = std::sqrt(m[2] * m[2] + 2 * m[2] * mass);
    m[1] = std::sqrt(m[1] * m[1] + 2 * m[1] * mass); 

    //    bull = new G4InuclElementaryParticle(bulletMomentum, bulletType); // counts mom[0] = E tot from mom[1]-mom[3]
    //    bull = new G4InuclParticle(bulletMomentum); // expects full mom[0]-mom[3] with correct E tot

    G4cout << G4endl << ">>> outgoing:" << G4endl;
    bull = new G4InuclElementaryParticle(m, 1); // expects full mom[0]-mom[3] with correct E tot
    bull->printParticle();

    G4cout << G4endl << ">>> Bullet initialization" << G4endl;

    m[3] = 0.0 ;
    m[2] = .585;
    m[1] = 0.0;
    m[0] = mass;

    G4InuclElementaryParticle* bull1 = new G4InuclElementaryParticle(m, 1); 
    bull1->printParticle();

    G4InuclParticle* bull2 = new G4InuclElementaryParticle(m, 1); 
    bull2->printParticle();
  }
  return 0;
}




  
