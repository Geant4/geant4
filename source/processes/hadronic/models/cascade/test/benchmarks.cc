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
// $Id: benchmarks.cc,v 1.22 2008-09-30 22:06:00 miheikki Exp $
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
#include "G4ios.hh"
#include "G4BertiniData.hh"
#include "G4CascadSpecialFunctions.hh"
#include "G4IonTable.hh"
#include "G4Nucleus.hh"
#include "G4NucleiModel.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "G4Proton.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4Fancy3DNucleus.hh"
#include "G4VParticleChange.hh"
#include "G4ParticleChange.hh"
#include "G4Track.hh"
#include "G4HadFinalState.hh"
#include "G4DynamicParticle.hh"

#include "G4CascadeInterface.hh"
#include "G4ElasticCascadeInterface.hh"
#include "G4InuclEvaporation.hh"

#include "G4CascadeKminusNChannel.hh"	   
#include "G4CascadeKminusPChannel.hh"
#include "G4CascadeKplusNChannel.hh"	   
#include "G4CascadeKplusPChannel.hh"
#include "G4CascadeKzeroBarNChannel.hh"	   
#include "G4CascadeKzeroBarPChannel.hh"
#include "G4CascadeKzeroNChannel.hh"	   
#include "G4CascadeKzeroPChannel.hh"


void test(std::string, int);
G4int tEvaporation(G4int A, G4int Z, G4double E);
G4int tCoulomb();
G4int benchmarkAll();
G4int tTiming();
G4int tBertiniData();
G4int tCrossSections();
G4int tInterface();
G4int tAnalyzer();
G4int tToyModel();
G4int tCascadeInterface();

int main(int argc, char **argv ) {
  G4int verboseLevel = 2;
  G4cout << "Geant4 cascade region benchmarks" << G4endl;

  if (argc < 2)
    {
      printf("usage: benchmarks <test ID> <parameters>\n");
      //      return(1);
    }
  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, foton = 10 };

  // defaults
  G4int testId      = 1;
  G4int nCollisions = 10;      // collisions to be generated
  G4int  bulletType = proton;    // bullet particle
  G4double     momZ = 160;      // momentum in z-direction
  G4double        A = 27.0;      // target atomic weight Al
  G4double        Z = 13.0;      // target atomic number
  testId      =           (argc > 2) ? atoi(argv[1]) : testId;
  nCollisions =           (argc > 2) ? atoi(argv[1]) : nCollisions;
  bulletType  =           (argc > 3) ? atoi(argv[2]) : proton;
  momZ        = G4double(((argc > 4) ? atoi(argv[3]) : momZ)) / GeV;
  A           = G4double(((argc > 5) ? atoi(argv[4]) : A));
  Z           = G4double(((argc > 7) ? atoi(argv[5]) : Z));

  if (verboseLevel > 1) {
    G4cout << " nCollisions " << nCollisions << G4endl;
    G4cout << "  bulletType " << bulletType  << G4endl;
    G4cout << "        momZ " << momZ        << G4endl;
    G4cout << "           A " << A           << G4endl;
    G4cout << "           Z " << Z           << G4endl;
  }

  //    test("Cascade interface", tCascadeInterface());

    // test("Evaporation", tEvaporation(27,13, 100));
  //test("Coulomb",               tCoulomb());

  //  test("Evaporation interface", tEvaporation());

  /*
    test("Timing",                tTiming());

    test("Bertini data",          tBertiniData());
    test("Cross sections",        tCrossSections());
    test("Toy model",             tToyModel());
    test("Cascade interface",     tCascadeInterface());
    test("Interface",             tInterface());
    test("Productions",           benchmarkAll()); // Run all models in tandem
    test("Analyzer",              tAnalyzer());  

  */
    test("Analyzer",              tAnalyzer());  
  return 0;       
}

void test(std::string txt, int testStatus) {
  G4cout << txt << ": ";
  if (testStatus){ 
    G4cout << "OK";
  } else {
    G4cout << "Fail" << G4endl;
  }; 

  G4cout << G4endl;  // test timing 
}


// Test program for G4 Bertini Evaporation.
G4int tElastic() { // test G4lasticCascadeInterface.hh
  G4cout << ">>> tElastic start" << G4endl;
  //  G4ElasticCascadeInterface *elastic = new G4ElasticCascadeInterface();
  G4cout << ">>> tElastic end" << G4endl;

  return 0;
}

int tEvaporation(G4int A, G4int Z, G4double E) { // test evaporation
  //G4LayeredNucleus nucl;
  G4InuclEvaporation bert;
  G4VParticleChange * pc;
  G4cout << ">>> tEvaporation start" << G4endl;
  G4cout << A << " " << Z << " " << E << G4endl;
  // nucl.SetParameters( A, Z);
  //nucl.AddExcitationEnergy( E );

  //  bert.setVerboseLevel(0);
  
  for ( G4double energy = 1 ; energy < 150 ; energy += 0.5 )
    {
         G4int n=0, p=0, d=0, t=0, h3=0, h4=0, g=0;

      //    nucl.SetParameters( A, Z);
      //nucl.AddExcitationEnergy( E - nucl.GetEnergyDeposit() );

      //pc = bert.BreakItUp( nucl );

      /*
	for ( G4int i = 0 ; i < pc->GetNumberOfSecondaries() ;  i++ )
	{
	char * name = pc->GetSecondary( i )->GetDefinition()->GetParticleName() ;
	if ( strcmp ( name , "proton" ) == 0 ) p++;
	if ( strcmp ( name , "neutron" ) == 0 ) n++;
	if ( strcmp ( name , "deuteron" ) == 0 ) d++;
	if ( strcmp ( name , "triton" ) == 0 ) t++;
	if ( strcmp ( name , "He3" ) == 0 ) h3++;
	if ( strcmp ( name , "alpha" ) == 0 ) h4++;
	if ( strcmp ( name , "gamma" ) == 0 ) g++;
	} // loop over particle change vector

	G4cout << energy << "\t"
	<< n << " "
	<< p << " "
	<< d << " "
	<< t << " "
	<< h3 << " "
	<< h4 << " "
	<< g << endl;
      */

    } // energy loop
  
  G4cout << ">>> tEvaporation end" << G4endl;
  return 1;
}

G4int tCoulomb()  {
  G4int verboseLevel = 1; // For benchmarking  quals 1.
  
  const G4int to_report = 1;
  G4int    nrain = 100;  // number of interactions to be generated
  G4double eMin  = 0.1;  // minimun energy for bullet
  G4double eMax  = 10.0; // maximum energy for bullet
  G4int    eBins = 10;   // bullet energy bins
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

    // proton momentum in Z-direction [GeV]
    G4double bulletEnergy = eMin + eStep * e;
 
    if (verboseLevel > 1) {
      G4cout << "Bullet E =" << bulletEnergy << " GeV" << G4endl;
    };

    G4InuclParticle* bull = new G4InuclElementaryParticle(bulletEnergy, 1); 
    G4InuclParticle* targ = NULL;

    switch (e) {
      // ::: add standard H,Be, Cu, Lb, U
    case 1:
      targ = new G4InuclNuclei(0.0, 197.0, 79.0);     // Au197 target at rest
    default:
      targ = new G4InuclNuclei(0.0, 208.0, 82.0); // Pb
    }

#ifdef CHECK_MOMC
    std::vector<G4double> total_mom_in = bull->getMomentum();
    std::vector<G4double> momt = targ->getMomentum();
    for(G4int i = 0; i < 4; i++) total_mom_in[i] += momt[i];
    std::vector<G4double> total_mom_out;
    bull->printParticle();
    targ->printParticle();
    if (verboseLevel > 1) {
      G4cout <<std::setw(15)<< " tot in mom: px " << total_mom_in[1] << " py " << total_mom_in[2] << " pz " << total_mom_in[3] << " e " << total_mom_in[0] << G4endl;
    }

#endif
    for(G4int i = 0; i < nrain; i++) {
      if((i + 1) % to_report == 0) 
	if (verboseLevel > 1) {
	  G4cout << " Event " << i+1 <<":" << G4endl;
      	}

      G4CollisionOutput cascadeParticles = collider->collide(bull, targ); // standard method

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
      analyser->printResults();
    //  analyser->printResultsSimple();
    //analyser->printResultsNtuple();
  }

  return 1;
}

G4int benchmarkAll()  {
  G4int verboseLevel = 1; // For benchmarking  quals 1.
  
  const G4int to_report = 1;
  G4int    nrain = 100;  // number of interactions to be generated
  G4double eMin  = 0.1;  // minimun energy for bullet
  G4double eMax  = 10.0; // maximum energy for bullet
  G4int    eBins = 10;   // bullet energy bins
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

    // proton momentum in Z-direction [GeV]
    G4double bulletEnergy = eMin + eStep * e;
 
    if (verboseLevel > 1) {
      G4cout << "Bullet E =" << bulletEnergy << " GeV" << G4endl;
    };

    G4InuclParticle* bull = new G4InuclElementaryParticle(bulletEnergy, 1); 
    G4InuclParticle* targ = NULL;


    switch (e) {
      // ::: add standard H,Be, Cu, Lb, U
    case 1:
      targ = new G4InuclNuclei(0.0, 197.0, 79.0);     // Au197 target at rest
    default:
      targ = new G4InuclNuclei(0.0, 208.0, 82.0); // Pb
    }

#ifdef CHECK_MOMC
    std::vector<G4double> total_mom_in = bull->getMomentum();
    std::vector<G4double> momt = targ->getMomentum();
    for(G4int i = 0; i < 4; i++) total_mom_in[i] += momt[i];
    std::vector<G4double> total_mom_out;
    bull->printParticle();
    targ->printParticle();
    if (verboseLevel > 1) {
      G4cout <<std::setw(15)<< " tot in mom: px " << total_mom_in[1] << " py " << total_mom_in[2] << " pz " << total_mom_in[3] << " e " << total_mom_in[0] << G4endl;
    }

#endif
    for(G4int i = 0; i < nrain; i++) {
      if((i + 1) % to_report == 0) 
	if (verboseLevel > 1) {
	  G4cout << " Event " << i+1 <<":" << G4endl;
      	}

      G4CollisionOutput cascadeParticles = collider->collide(bull, targ); // standard method

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

  return 1;

}

int tTiming() { // Test speed of pow(x, 2) compared to x * x
  G4double y[3]= {1.0, 1.3, 1.2}; 
  y[1]=1.2;
  G4int LOOPS = 20000000; // Set test parameters

  clock_t startTime;
  clock_t endTime;

  startTime = clock();
  G4double ans;
  for(G4int i = 1; i < LOOPS; i++){
    ans = std::pow(y[2], 2);
  };
  G4cout << "ans =" << ans << G4endl;

  endTime = clock();
  //  G4double firstTime = (G4double)(endTime - startTime) /
  //  (CLOCKS_PER_SEC * 1000000.0);

  G4double firstTime = (G4double)(endTime - startTime);

  G4cout << "pow(x, 2) time: " << firstTime  << G4endl;

  startTime = clock();
  for(G4int j = 1; j < LOOPS; j++){
    ans = y[2] * y[2];
  };
  G4cout << "ans =" << ans << G4endl;
  endTime = clock();
  //  G4double secondTime = (G4double)(endTime - startTime) / 
  //  (CLOCKS_PER_SEC * 1000000.0);


  G4double secondTime = (G4double)(endTime - startTime);
  G4cout << "x * x time: " << secondTime << G4endl;
  G4cout << "pow / * speed ratio = " << firstTime / secondTime << G4endl;
  return 1;
}



int tBertiniData() // test and demontrate singleton usage
{
  G4cout << G4endl << "testing G4BertiniData" << G4endl;
  G4BertiniData *db = new G4BertiniData(); // sever

  // G4BertiniData *data1 = db->Instance(); // client 
  //G4BertiniData *data2 = db->Instance(); // old instance used
  delete db;
  return 1;
}    


G4int tCrossSections() {   // print cross section data
  G4int verboseLevel = 1;

  if (verboseLevel > 1) {
    G4cout << " >>> crossSections() " << G4endl;
  }

  if (verboseLevel > 1) {
    G4cout << " MeV: " << MeV << " GeV: " << GeV << G4endl;
  }

  // 100  <> 10 GeV
  const G4int types[] = { 1, 2, 3, 5, 7};

  for (G4int iE = 1; iE < 151; iE++) { 

    if (verboseLevel > 0) {
      G4cout.precision(4);
      G4double e = G4double(iE) / 10.0;
      G4cout << std::setw(9)  << e;

      for (G4int j = 0; j < 5; j++) {
	G4cout << std::setw(9)  << G4CascadSpecialFunctions::crossSection(e, types[j]);
      }}
    G4cout << G4endl;
  }
  return 1;
}

int tInterface() { // test iterface

  typedef std::vector<G4InuclElementaryParticle>::iterator particleIterator;
  typedef std::vector<G4InuclNuclei>::iterator nucleiIterator;

  enum particleType { nuclei = 0, proton = 1, neutron = 2, pionPlus = 3, pionMinus = 5, pionZero = 7, photon = 10 };
  G4int verboseLevel = 3;

  if (verboseLevel > 1) {
    G4cout << "tInterface::   begin" << G4endl; 
  }

  //  G4int bulletType = 0;

  std::vector<G4double>  momentumBullet(4, 0.0);
  momentumBullet[0] = 1.37126;
  momentumBullet[3] = 1.5;

  G4InuclParticle *  bullet = new G4InuclElementaryParticle(momentumBullet, 1); 

  G4InuclNuclei*   target  = NULL;
  G4InuclParticle* targetH = NULL;
  G4NucleiModel*   model = NULL;

  G4double theNucleusA = 1;
  std::vector<G4double> targetMomentum(4, 0.0);

  G4CollisionOutput output;

  G4ElementaryParticleCollider*   colep = new G4ElementaryParticleCollider;
  G4IntraNucleiCascader*            inc = new G4IntraNucleiCascader; // the actual cascade
  inc->setInteractionCase(1); // Interaction type is particle with nuclei.

  G4NonEquilibriumEvaporator*     noneq = new G4NonEquilibriumEvaporator;
  G4EquilibriumEvaporator*         eqil = new G4EquilibriumEvaporator;
  G4Fissioner*                     fiss = new G4Fissioner;
  G4BigBanger*                     bigb = new G4BigBanger;

  G4InuclCollider*             collider = new G4InuclCollider(colep, inc, noneq, eqil, fiss, bigb);

  if (verboseLevel > 1) {
    G4cout << "tInterface::   Bertini classes instantiated" << G4endl; 
  }

  for (G4int i = 1; i< 100 ; i++) {
    if ( theNucleusA < 1.5 ) 
      {
	model = new G4NucleiModel(new G4InuclNuclei(targetMomentum, 1, 1));
	targetH = new G4InuclElementaryParticle((model->generateNucleon(1, 1)).getMomentum(), 1); 
   
	//		do
	{
	  G4cout << "+";
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

    if (verboseLevel > 1) {
      G4cout << "tInterface::   collide done" << G4endl; 
    }
    std::vector<G4InuclNuclei>             nucleiFragments = output.getNucleiFragments();
    std::vector<G4InuclElementaryParticle> particles =       output.getOutgoingParticles();

    G4int numSecondaries = nucleiFragments.size()+particles.size();
    G4cout << "num secondaries: " << numSecondaries << G4endl;
    if(!particles.empty()) { 
      particleIterator ipart;
      G4int outgoingParticle;

      G4double ekin;
      for(ipart = particles.begin(); ipart != particles.end(); ipart++) {
	outgoingParticle = ipart->type();
	std::vector<G4double> mom = ipart->getMomentum();
	ekin = ipart->getKineticEnergy() * GeV;
	G4ThreeVector aMom(mom[1], mom[2], mom[3]);
	aMom = aMom.unit();

	if (verboseLevel > 2) {
	  G4cout << "ekin" << ekin << G4endl; 
	}
      
      }

      nucleiIterator ifrag;

      for(ifrag = nucleiFragments.begin(); ifrag != nucleiFragments.end(); ifrag++) 
	{
	  ekin = ifrag->getKineticEnergy() * GeV;
	  std::vector<G4double> mom = ifrag->getMomentum();
	  G4ThreeVector aMom(mom[1], mom[2], mom[3]);
	  aMom = aMom.unit();

	  // hpw @@@ ==> Should be zero: G4double fragmentExitation = ifrag->getExitationEnergyInGeV();


	  G4int A = G4int(ifrag->getA());
	  G4int Z = G4int(ifrag->getZ());

	  if (verboseLevel > 2) {
	    G4cout << " Nuclei fragment: " << G4endl;
	    G4cout << "A: " << A << " Z: " << Z << G4endl;
	    ifrag->printParticle();
	  }
	}
    }
  }

  if (verboseLevel > 1) {
    G4cout << "tInterface::   end" << G4endl; 
  }
  return 1;
}

int tAnalyzer() { 
  //  G4int verboseLevel = 1;
  // General test program for hetc and inucl

  const G4int to_report = 1;
  G4int nrain = 1000; // number of interactions to be generated
  G4double eMin  = 1.0; // minimun energy for bullet
  G4double eMax  = 2.0;   // maximum energy for bullet
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
    G4cout << "Bullet E =" << bulletEnergy << " GeV" << G4endl;
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
    std::vector<G4double> total_mom_in = bull->getMomentum();
    std::vector<G4double> momt = targ->getMomentum();
    for(G4int i = 0; i < 4; i++) total_mom_in[i] += momt[i];
    std::vector<G4double> total_mom_out;
    bull->printParticle();
    targ->printParticle();
    G4cout << " tot in mom: px " << total_mom_in[1] << " py " << total_mom_in[2]
	   << " pz " << total_mom_in[3] << " e " << total_mom_in[0] << G4endl;
#endif
    for(G4int i = 0; i < nrain; i++) {
      if((i + 1) % to_report == 0) 
	//	if (verboseLevel > 1) {
	G4cout << " Event " << i+1 <<":" << G4endl;
      //	}
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
      analyser->printResults();
    //  analyser->printResultsSimple();
    //    analyser->printResultsNtuple();
  }

  return 1;       
}


int tToyModel() {
  G4cout << G4endl << "       ... intra-nuclear cascade toy model ..." << G4endl;
	 
  // Model based on: 
  // Methods in Subnuclear Physics, Vol IV, Part 3 pp. 137-144,
  // Proceedings of the International School of Elelmantary Particle Physics, 
  // Herceg-Novi, Yugoslavia, 1968   
  // Ed. M. Mikolic

  // Set constanst and create terminology
  G4double cutOff = 20 * MeV;
  G4double fm = 1.0 * std::pow(10.0, -15);
  G4double oneThird = 0.3333333333333;
  enum particleType {proton, neutron};
  // Eenum channelType {pionProduction, };
  enum collisionType {elastic, inelastic};

  std::vector<G4KineticTrack*>  particleVector;   

  // Set test parameters
  G4int        numberOfCascades       = 3;
  G4int        targetNucleusN         = 10;
  G4int        targetNucleusZ         = 10;
  particleType incidentParticleType   = proton;
  G4double     incidentParticleEnergy = 1000.0 * MeV;                         
  G4Fancy3DNucleus      theCarbon;
  G4double              theA = 12.0;
  G4double              theZ =  6.0;
  theCarbon.Init(theA, theZ);

  G4cout << G4endl << " Nucleus mass number = " << theCarbon.GetMassNumber() << G4endl
	 << " Nucleus mass = " << theCarbon.GetMass() << G4endl
	 << " Nucleus charge = " << theCarbon.GetCharge() << G4endl;
 
  G4cout << " projectile                 : " << incidentParticleType             << G4endl;
  G4cout << " projectile energy          : " << incidentParticleEnergy << " MeV" << G4endl;
  G4cout << " neutrons in target nucleus : " << targetNucleusN                   << G4endl;
  G4cout << " protons  in target nucleus : " << targetNucleusZ                   << G4endl;

  for (G4int cascadeIndex = 1; cascadeIndex <= numberOfCascades; cascadeIndex++){

    G4cout << G4endl << "........................ cascade " << cascadeIndex << 
      " .............................." << G4endl;

    // Select the impact parameter b
    G4double targetNucleusA = targetNucleusN + targetNucleusZ;
    G4double radius0        = 1.0 * fm;                                       // radius of hydrogen atom nucleus
    G4double bMax           = radius0 * fm / std::pow(targetNucleusA, oneThird) ;  // radius of atom with a = n + z
    G4double b              = bMax * std::sqrt(G4UniformRand());
    G4double z              = b;                                              // x, y, z coordinates of impact point
    G4double x              = - std::sqrt(sqr(bMax) - sqr(z));
    G4double y              = 0;
    //    G4cout << radius0 << " " << bMax << endl;
    G4cout << "impact point (x, y, z) : "  << "\t" << x <<", " << "\t" << y << ", " << "\t" << z  << G4endl;

    // *) Find the interaction distance b
    G4double xSecNeutron = 1 * millibarn;
    G4double xSecProton  = 1 * millibarn;
    G4double lambda      = 4 * pi * targetNucleusA * std::pow(radius0, 3) /
      (3 * (targetNucleusZ * xSecProton + targetNucleusN * xSecNeutron));
    b = - lambda * std::log(G4UniformRand());

    G4cout << "interaction distance   : " << "\t" << b << G4endl;

    // Choose the sruck particle (n/p)
    particleType targetParticle;
    if (G4UniformRand() < targetNucleusZ * xSecProton / 
	( targetNucleusZ * xSecProton + targetNucleusN * xSecNeutron))
      targetParticle = proton;
    else 
      targetParticle = neutron;

    G4cout << "target particle        : " << "\t" << targetParticle << G4endl;

    // Choose the struck nucleon momentum
    //    G4double qMax = (300.0 / radius0) * std::pow(10, -15) * MeV;
    G4double qMax = (300.0 / radius0) * pow(10, -15) * MeV;
    G4double q    = qMax * std::sqrt(G4UniformRand());

    G4cout << "target  momentum       : " << "\t" << q << G4endl;

    // Choose the direction of the struck nucleon
    G4double mu  = 1 - 2 * G4UniformRand();
    G4double phi = 2 * pi * G4UniformRand();
    G4cout << "target direction       : " << "\t" << mu << ", " << "\t" << phi << G4endl;

    // Choose the nature of the collision (elastic/inelastic)
    collisionType collision;
    if (G4UniformRand() < 0.8) // replace with more realistic 
      collision = elastic;
    else
      collision = inelastic;
    G4cout << "collision type         : " << "\t" << collision << G4endl;

    if (collision == elastic) {
      // Choose the scattering angle theta
      G4double thetaCos = 2 *G4UniformRand() - 1;                           // dummy
      G4cout << " scattering angle      : " << "\t" << thetaCos << G4endl;
    }

    if (collision == inelastic) {

      G4ParticleDefinition* theProtonDefinition = G4Proton::ProtonDefinition();
      G4double              theProtonFormationTime = 0.0*ns;
      G4ThreeVector         theProtonInitialPosition(0.0 * cm, 0.0 * cm, 0.0 * cm);
      G4LorentzVector       theProton4Momentum(0.0 * MeV, 0.0 * MeV, 1.0 * MeV, 938.27 * MeV);
      G4KineticTrack        aProton(theProtonDefinition, 
				    theProtonFormationTime,
				    theProtonInitialPosition,
				    theProton4Momentum);

      G4KineticTrack     *xx = new G4KineticTrack(theProtonDefinition, 
						  theProtonFormationTime,
						  theProtonInitialPosition,
						  theProton4Momentum);
                          
      G4Fancy3DNucleus      theCarbon;
      G4double              theA = 12.0;
      G4double              theZ =  6.0;
      theCarbon.Init(theA, theZ);
   
      particleVector.push_back(xx);
      G4cout << " pion production       : " << G4endl;

    }

    // Choose the azimuthal angle of rotation
    G4double alpha = 2 * pi * G4UniformRand();
    G4cout << " azimuthal angle       : " << "\t" << alpha << G4endl;

    // Transform the energy back to the laboratory system (lab)
    G4double energy =  50 * MeV* G4UniformRand();                            //dummy

    // For all items in particleVectos
    // Pauli exclusion principle
    G4double M = 0.5 ;                                                       //dummy
    G4double fermiEnergy = sqr(qMax)/(2 * M);
    if (energy <  fermiEnergy)  {                                            // state occupied
      G4cout << "pauli blocking " << G4endl;
      // Add article to traced list  
    }          

    // Transform the allowed collision products to the laboratory system

    // Follow the products until the escape the nucleus or energy is below cut-off
    if (energy < cutOff) {
      G4cout << "energy below cut-off   "<< G4endl;
      // Add energy to atom exitation
      // exitation += particle.energy 
      // update atom quantum numbers
      // if particle = proton, a += 1, z += 1 
      // if particle = neutron, a += 1, n += 1 
    }
  }
   
  return 1;
}    
/*
  ... intra-nuclear cascade toy model ...

  Nucleus mass number = 12
  Nucleus mass = 11174.9
  Nucleus charge = 6
  projectile                 : 0
  projectile energy          : 1000 MeV
  neutrons in target nucleus : 10
  protons  in target nucleus : 10

  ........................ cascade 1 ..............................
  impact point (x, y, z) :        -2.82281e-31,   0,      2.36724e-31
  interaction distance   :        3.54704e-20
  target particle        :        1
  target  momentum       :        152.393
  target direction       :        -0.713471,      2.86607
  collision type         :        0
  scattering angle      :        0.159015
  azimuthal angle       :        4.14197
  pauli plocking

  ........................ cascade 2 ..............................
  impact point (x, y, z) :        -2.82281e-31,   0,      2.36724e-31
  interaction distance   :        3.54704e-20
  target particle        :        1
  target  momentum       :        152.393
  target direction       :        -0.713471,      2.86607
  collision type         :        0
  scattering angle      :        0.159015
  azimuthal angle       :        4.14197
  pauli plocking
  energy below cut-off  

*/

int tCascadeInterface() {
  G4int verboseLevel = 2;                          
  G4int numberOfCascades = 2;
 
  if (verboseLevel > 1) G4cout << ">>> tCascadeInterface start" << G4endl;

  //G4ParticleDefinition *particle = G4PionMinus::PionMinus();  
  //G4ParticleDefinition *particle = G4Neutron::Neutron();  
  G4ParticleDefinition *particle = G4Proton::Proton();  
  G4DynamicParticle *projectile = new G4DynamicParticle(); 
  projectile->SetDefinition(particle);
  //projectile->SetKineticEnergy( 1.0 * GeV);
  G4double projectileMomentum = 1.0 * GeV;
  projectile->SetMomentum(projectileMomentum);
  projectile->SetMomentumDirection(1.0, 0.0, 0.0);  

  if (verboseLevel > 1) {
    G4cout << "projectile" << G4endl;
    G4cout << " type           : " << projectile->GetDefinition()        << G4endl;
    G4cout << " kinetic energy : " << projectile->GetKineticEnergy()     << G4endl;
    G4cout << " momentum       : " << projectile->GetMomentum()          << G4endl;
    G4cout << " p direction    : " << projectile->GetMomentumDirection() << G4endl;
  }
  G4Nucleus targetNucleus;                                        
  G4double a(10);
  G4double z(10);
  targetNucleus.SetParameters(a, z);

  if (verboseLevel > 1) {
    G4cout << "target" << G4endl;
    G4cout << " a              : " << a                              << G4endl;
    G4cout << " z              : " << z                              << G4endl;
    G4cout << " atomic mass    : " << targetNucleus.AtomicMass(a, z) << G4endl;
  }
  G4ThreeVector       inVector(0.0, 0.0, 2000.0);
  G4ThreeVector       outVector;
  G4Proton          * aProton   = G4Proton::Proton();
  G4DynamicParticle   aParticle;
  aParticle.SetDefinition(aProton);
  G4double Momentum = 2000;
  G4double Tkin = sqrt(Momentum*Momentum+938.27*938.27)-938.27;
  inVector.setZ(Momentum);
  G4cout <<  "  inVector " << inVector.x() << " " << inVector.y() << " " << inVector.z() << " " << " Tkin " << Tkin << G4endl;
  aParticle.SetMomentum(inVector);
  aParticle.SetKineticEnergy(Tkin);
  const  G4HadProjectile hadProj = G4HadProjectile(aParticle);
  G4HadFinalState * hadSta = new G4HadFinalState();   
  G4CascadeInterface *theCascade  = new G4CascadeInterface();  // instantiate Bertini cascade interface
  for (G4int cascadeID = 1; cascadeID <= numberOfCascades; cascadeID++) { 
    if (verboseLevel > 1) G4cout << "inc " << cascadeID << G4endl;
    hadSta = theCascade->ApplyYourself(hadProj, targetNucleus); // apply interface
    G4int nPart = hadSta->GetNumberOfSecondaries(); G4cout << "  # secondaries " << nPart << G4endl;
    outVector = hadSta->GetMomentumChange(); G4cout << "  momentum change " << outVector << G4endl;
    G4double outE = hadSta->GetEnergyChange(); G4cout << "  energy change " << outE << G4endl;
    
    for (G4int iSecondary =1 ; iSecondary < nPart; iSecondary++) { 
      G4cout << " secondary " << iSecondary << G4endl;  
      G4HadSecondary * NuclSecond = hadSta->GetSecondary(iSecondary);
      G4DynamicParticle * secPart = NuclSecond->GetParticle(); G4cout<<" nucleus name " << secPart->GetDefinition()->GetParticleName() << G4endl;
      G4ThreeVector outVectorN =  secPart->GetMomentum(); G4cout << " out vector " << outVectorN << G4endl;
      G4double outEtot = secPart->GetTotalEnergy(); G4cout << " particle  tot E " << outEtot << G4endl;  
      G4double outEkin = secPart->GetKineticEnergy(); G4cout << " particle  kin E " << outEkin << G4endl;  
    }
  }
  delete projectile;
  delete theCascade;
  return 1;   
}
