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
// $Id: evaporation.cc,v 1.1 2007-11-13 23:45:43 miheikki Exp $
// based on V.Lara evaporation test
#include <iostream>
#include <iomanip>
#include "globals.hh"
#include "Randomize.hh"

#include "G4InuclEvaporation.hh" // evaporations to be tested
#include "G4Evaporation.hh"

#include "G4NucleiProperties.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Electron.hh"
#include "G4Fragment.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ThreeVector.hh" 

#include <vector>

#include <iomanip>
#include <iostream>
#include <string>
#include <stdio.h>
#include <time.h>

//
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


int main() {

  /*
    G4ParticleDefinition *theGamma = G4Gamma::GammaDefinition();
    G4ParticleDefinition *theElectron = G4Electron::ElectronDefinition();
    G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
    G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();   
    G4ParticleDefinition *theDeuteron = G4Deuteron::DeuteronDefinition();
    G4ParticleDefinition *theTriton = G4Triton::TritonDefinition();
    G4ParticleDefinition *theHelium3 = G4He3::He3Definition();
    G4ParticleDefinition *theAlpha = G4Alpha::AlphaDefinition();
  */

  G4cout << G4endl << "Evaporation Test program" << G4endl;
G4cout << "Test Momentum  four vectors should be zero.[MeV]" << G4endl; 
  G4InuclEvaporation * theEvaporation = new G4InuclEvaporation;
  //  theEvaporation->setVerboseLevel(10);
  //  G4Evaporation * theEvaporation = new G4Evaporation;

  G4int mode = 1;

  while (mode != 0 && mode != 1) {
    G4cout << "Mode (0 or 1): ";
    G4cin >> mode;
  }



  if (mode == 0) {
   
    G4cout << "A = ";
    G4int MyA;
    G4cin >> MyA;
    G4cout << "Z = ";
    G4int MyZ;
    G4cin >> MyZ;
    G4cout << "ExcitationEnergy (MeV) = ";
    G4double MyExE;
    G4cin >> MyExE;
    //    MyExE *= MyA;
    G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV + MyExE;
    G4cout << "Momentum (MeV): " << G4endl;
    G4cout << "                Px = ";
    G4double MyPx;
    G4cin >> MyPx;
    G4cout << "                Py = ";
    G4double MyPy;
    G4cin >> MyPy;
    G4cout << "                Pz = ";
    G4double MyPz;
    G4cin >> MyPz;
    
    G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
    //    G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
    G4LorentzVector initialMomentum(triV,std::sqrt(triV.mag2()+AtomicMass*AtomicMass*MeV*MeV));

    // put info about excited nucleus in fragment class
    G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);

    G4int events = 1;
    G4cout << "Iterations: ";
    G4cin >> events;

    G4double NofP = 0.0;
    G4double NofN = 0.0;

    for (G4int i = 0; i < events; i++) {
      G4cout << "Iteration: " << i+1 << G4endl;
      G4cout << "----------------" << G4endl;
      G4cout << "     Initial fragment" << G4endl;
      // G4cout << theExcitedNucleus << G4endl << G4endl;

      G4cout << "     Fragments evaporated" << G4endl << G4endl;
      // DeExcite the nucleus 
      G4FragmentVector * theFragVector = theEvaporation->BreakItUp(theExcitedNucleus);

      //G4cout << "#fragments " << theFragVector->size()<< G4endl;

      G4LorentzVector TestMomentum(initialMomentum);
      for (G4int j=0; j < (int)theFragVector->size(); j++) {
	//	G4cout << theFragVector->at(j) << G4endl;

	// Test 4-momentum conservation 
	TestMomentum -= theFragVector->at(j)->GetMomentum();

	// Calculating multiplicities
	if (theFragVector->at(j)->GetA() == 1 && theFragVector->at(j)->GetZ() == 1) NofP++;
	else if (theFragVector->at(j)->GetA() == 1 && theFragVector->at(j)->GetZ() == 0) NofN++;
      }

      G4cout << "******************" << G4endl;
      G4cout << "* Test Momentum = " << TestMomentum << G4endl;
      G4cout << "******************" << G4endl;
      
      theFragVector->clear();
      delete theFragVector;
    }

    G4cout << "Multiplicities: Neutrons -> " << NofN/events << " Protons -> " << NofP/events << G4endl;
  } else if (mode == 1) {


    G4int events = 100;
    G4cout << "Number of events: " << events << G4endl;
 
    for (G4int i = 0; i < events; i++) {

      //            G4cout << "Event number: " << i+1 << G4endl;
      G4cout << "Ev" << i+1 <<  " " ;
      //G4cout << "--------------------" << G4endl;
      
      /*
	G4int MyA = RandFlat::shoot(17,200);
	G4int MyZ = RandFlat::shoot(G4NucleiPropertiesTable::MinZ(MyA), G4NucleiPropertiesTable::MaxZ(MyA));
	G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV;      
	G4double MyExE = RandFlat::shoot(2.0*MyA,10.0*MyA);
	G4double MyPx = RandFlat::shoot(-2000,2000);
	G4double MyPy = RandFlat::shoot(-2000,2000);
	G4double MyPz = RandFlat::shoot(-2000,2000);
      */
      G4int MyA = 197;
      G4int MyZ = 79;

      G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV;      
      //      G4double MyExE = 10+G4UniformRand()*1000;
      G4double MyExE = 3000;
      G4double MyPx = 10;
      G4double MyPy = 100;
      G4double MyPz = 1000;

      G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
      G4LorentzVector initialMomentum(triV,std::sqrt(triV.mag2()+
						     (AtomicMass*MeV+MyExE*MeV)*
						     (AtomicMass*MeV+MyExE*MeV))
				      );
    
      //      G4cout << ":::::MeV:e1" << MeV;
      // put info about excited nucleus in fragment class
      G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);
      
      //      G4cout << "Excited fragment: "<< G4endl;
      //G4cout << theExcitedNucleus << G4endl;

      G4FragmentVector * theFragVector = theEvaporation->BreakItUp(theExcitedNucleus);

      //G4cout << "test: #fragments: " << theFragVector->size()<< G4endl;


       G4LorentzVector TestMomentum(initialMomentum);
      for (G4int j=0; j < (int)theFragVector->size(); j++) {
	//	G4cout << ":::::::::::" << theFragVector->at(j)->GetMomentum() << G4endl;
       	TestMomentum -= theFragVector->at(j)->GetMomentum(); 	// Test 4-momentum conservation 
      }

      G4ThreeVector t(0,0,0);
    //    G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
    G4LorentzVector initialM(t, MyExE);


    G4cout << TestMomentum - initialM << " [MeV]" << G4endl; // should be (0,0,0,0)
    // G4cout << "******************" << G4endl;
    //    G4cout << "* Test Momentum = " << TestMomentum - initialM << " [MeV]" << G4endl; // should be (0,0,0
    //  G4cout << "******************" << G4endl;
      
      //      G4cout << "exitation  was E: " << MyExE << G4endl;
      theFragVector->clear();
      delete theFragVector; 
    }  
  }

  delete theEvaporation;
  return 0;
}
