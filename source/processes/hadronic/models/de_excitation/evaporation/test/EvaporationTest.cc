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
//

#include <iostream>
#include <iomanip>
#include "globals.hh"

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


int main()
{
  G4ParticleDefinition *theGamma = G4Gamma::GammaDefinition();
  G4ParticleDefinition *theElectron = G4Electron::ElectronDefinition();
  G4ParticleDefinition *theNeutron = G4Neutron::NeutronDefinition();
  G4ParticleDefinition *theProton = G4Proton::ProtonDefinition();   
  G4ParticleDefinition *theDeuteron = G4Deuteron::DeuteronDefinition();
  G4ParticleDefinition *theTriton = G4Triton::TritonDefinition();
  G4ParticleDefinition *theHelium3 = G4He3::He3Definition();
  G4ParticleDefinition *theAlpha = G4Alpha::AlphaDefinition();


  theProton->SetCuts(1.0);
  theGamma->SetCuts(1.0);
  theElectron->SetCuts(1.0);
  theNeutron->SetCuts(1.0);
  theDeuteron->SetCuts(1.0);
  theTriton->SetCuts(1.0);
  theHelium3->SetCuts(1.0);
  theAlpha->SetCuts(1.0);


  cout << G4endl << "Evaporation Test program" << G4endl;
 
  G4Evaporation * theEvaporation = new G4Evaporation;

  G4int mode = -1;
  while (mode != 0 && mode != 1) {
    cout << "Mode (0 or 1): ";
    G4cin >> mode;
  }

  cout << G4endl << G4endl << G4endl;

  if (mode == 0) {
   
    cout << "A = ";
    G4int MyA;
    G4cin >> MyA;
    cout << "Z = ";
    G4int MyZ;
    G4cin >> MyZ;
    cout << "ExcitationEnergy (MeV) = ";
    G4double MyExE;
    G4cin >> MyExE;
    //    MyExE *= MyA;
    G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV + MyExE;
    cout << "Momentum (MeV): " << G4endl;
    cout << "                Px = ";
    G4double MyPx;
    G4cin >> MyPx;
    cout << "                Py = ";
    G4double MyPy;
    G4cin >> MyPy;
    cout << "                Pz = ";
    G4double MyPz;
    G4cin >> MyPz;
    

    G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
    //    G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
    G4LorentzVector initialMomentum(triV,std::sqrt(triV.mag2()+AtomicMass*AtomicMass*MeV*MeV));
    

    // put info about excited nucleus in fragment class
    G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);


    G4int events = 1;
    cout << "Iterations: ";
    G4cin >> events;

    G4double NofP = 0.0;
    G4double NofN = 0.0;

    for (G4int i = 0; i < events; i++) {

      cout << "Iteration: " << i+1 << G4endl;
      cout << "----------------" << G4endl;
      cout << "     Initial fragment" << G4endl;
      cout << theExcitedNucleus << G4endl << G4endl;

      cout << "     Fragments evaporated" << G4endl << G4endl;
      // DeExcite the nucleus 
      G4FragmentVector * theFragVector = theEvaporation->BreakItUp(theExcitedNucleus);

      
      G4LorentzVector TestMomentum(initialMomentum);
      for (G4int j=0; j < theFragVector->entries(); j++) {
	cout << theFragVector->at(j) << G4endl;

	// Test 4-momentum conservation 
	TestMomentum -= theFragVector->at(j)->GetMomentum();

	// Calculating multiplicities
	if (theFragVector->at(j)->GetA() == 1 && theFragVector->at(j)->GetZ() == 1) NofP++;
	else if (theFragVector->at(j)->GetA() == 1 && theFragVector->at(j)->GetZ() == 0) NofN++;
      }

      cout << "******************" << G4endl;
      cout << "* Test Momentum = " << TestMomentum << G4endl;
      cout << "******************" << G4endl;
      
      theFragVector->clearAndDestroy();
      delete theFragVector;
    }

    G4cout << "Multiplicities: Neutrons -> " << NofN/events << " Protons -> " << NofP/events << G4endl;
  } else if (mode == 1) {


    G4int events = 0;
    cout << "Number of events: ";
    G4cin >> events;

    for (G4int i = 0; i < events; i++) {

      cout << "Event number: " << i+1 << G4endl;
      cout << "--------------------" << G4endl;
      
      G4int MyA = RandFlat::shoot(17,200);

      G4int MyZ = RandFlat::shoot(G4NucleiPropertiesTable::MinZ(MyA),
				  G4NucleiPropertiesTable::MaxZ(MyA));

      G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV;
      
      G4double MyExE = RandFlat::shoot(2.0*MyA,10.0*MyA);
      G4double MyPx = RandFlat::shoot(-2000,2000);
      G4double MyPy = RandFlat::shoot(-2000,2000);
      G4double MyPz = RandFlat::shoot(-2000,2000);

      G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
      G4LorentzVector initialMomentum(triV,std::sqrt(triV.mag2()+
						(AtomicMass*MeV+MyExE*MeV)*
						(AtomicMass*MeV+MyExE*MeV))
				      );

    
      // put info about excited nucleus in fragment class
      G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);
      
      cout << "Excited fragment: "<< G4endl;
      cout << theExcitedNucleus << G4endl;
      
      G4FragmentVector * theFragVector = theEvaporation->BreakItUp(theExcitedNucleus);
      

      G4LorentzVector TestMomentum(initialMomentum);
      for (G4int j=0; j < theFragVector->entries(); j++) {
	//	cout << theFragVector->at(j) << G4endl;

	// Test 4-momentum conservation 
	TestMomentum -= theFragVector->at(j)->GetMomentum();
      }

      cout << "******************" << G4endl;
      cout << "* Test Momentum = " << TestMomentum << G4endl;
      cout << "******************" << G4endl;



      theFragVector->clearAndDestroy();
      delete theFragVector;
      
    }
    
  }

  delete theEvaporation;


  return 0;

}




