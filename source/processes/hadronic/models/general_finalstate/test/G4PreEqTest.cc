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
// * authors in the GEANT4 collaboration.                             *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
#include <iostream>
#include <iomanip>
#include "globals.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"


#include "G4PreCompoundModel.hh"
#include "G4ExcitationHandler.hh"
#include "G4NucleiProperties.hh"
#include "G4Fragment.hh"
#include "G4ThreeVector.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Electron.hh"
#include "G4GenericIon.hh"

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
G4ParticleDefinition *theIon =
		  G4GenericIon::GenericIonDefinition();
  theProton->SetCuts(1.0);
  theGamma->SetCuts(1.0);
  theElectron->SetCuts(1.0);
  theNeutron->SetCuts(1.0);
  theDeuteron->SetCuts(1.0);
  theTriton->SetCuts(1.0);
  theHelium3->SetCuts(1.0);
  theAlpha->SetCuts(1.0);


  G4cout << G4endl << "Pre-Equilibrium Test program" << G4endl;

  // HBOOK initialization
    HepTupleManager* HbookManager;
  HbookManager = new HBookFile("PreEqTest.hbook", 58);


  // Book a ntuple
  HepTuple* ntupleP;
  HepTuple* ntupleN;
//   ntuple = HbookManager->ntuple("Nucleon energy");

//   HepTuple* ntuple2;
//   ntuple2 = HbookManager->ntuple("NucleonMultiplicities");



  // Book histograms
//   HepHistogram* Protons;
//   Protons = HbookManager->histogram("Proton Mass", 100,930.,950.);

//   HepHistogram* Neutrons;
//   Neutrons = HbookManager->histogram("Neutron Mass", 100,930.,950.);


  



  G4ExcitationHandler * theExcitationHandlerPtr = new G4ExcitationHandler;
  
  G4PreCompoundModel * thePreCompoundModelPtr = new G4PreCompoundModel(theExcitationHandlerPtr);


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

    cout << "Number of excitation energies: ";
    G4int NumOfE;
    G4cin >> NumOfE;
    
    G4double * ExcitE = new G4double[NumOfE];
    for (G4int i = 0; i < NumOfE; i++) {
      cout << "Ecitation Energy " << i+1 << ": ";
      G4cin >> ExcitE[i];
    }

    //    cout << "ExcitationEnergy (MeV) = ";
//     G4double MyExE;
//     G4cin >> MyExE;
//     MyExE += G4NucleiPropertiesTable::GetBindingEnergy(MyZ,MyA)/(MyA*MeV);
//     //    MyExE *= MyA;
//     G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ)/MeV + MyExE;
    
    
    
   cout << "Momentum (MeV): " << G4endl;
    cout << "                Px = ";
    G4double MyPx = 0.0;
    G4cin >> MyPx;
    cout << "                Py = ";
    G4double MyPy = 0.0;
    G4cin >> MyPy; 
    cout << "                Pz = ";
    G4double MyPz = 0.0;
    G4cin >> MyPz;
    

    

    cout << "Number of Particles = ";
    G4int MyParticles;
    G4cin >> MyParticles;
    cout << "Number of Holes     = ";
    G4int MyHoles;
    G4cin >> MyHoles;
    cout << "Number of Charged   = ";
    G4int MyCharged;
    G4cin >> MyCharged;

    G4int events = 1;
    cout << "Iterations: ";
    G4cin >> events;

    
    
    for (G4int k=0; k < NumOfE; k++) {
      G4ThreeVector incidentProton(0.0,0.0,1.0);    
      G4double MyExE = ExcitE[k];
      
 //     MyPz = sqrt(MyExE*(MyExE + 2.0*theProton->GetPDGMass()/MeV));
      
      MyExE += G4NucleiPropertiesTable::GetBindingEnergy(MyZ,MyA)/(MyA*MeV); 
      G4double AtomicMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(MyZ,MyA)/MeV + MyExE;

      ntupleP = HbookManager->ntuple("Protons");
      ntupleN = HbookManager->ntuple("Neutrons");

      G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
      //    G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
      G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass*MeV*MeV));
    

      // put info about excited nucleus in fragment class
      G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);
      theExcitedNucleus.SetNumberOfExcitons(MyParticles+MyHoles);
      theExcitedNucleus.SetNumberOfHoles(MyHoles);
      theExcitedNucleus.SetNumberOfCharged(MyCharged); 

//       G4cout << "The initial Fragment: " << theExcitedNucleus << G4endl;
      
      G4double SumNP = 0;
      G4double SumNN = 0; 
      
      for (G4int i = 0; i < events; i++) {
	
	//	cout << "Iteration: " << i+1 << G4endl;
	//	cout << "----------------" << G4endl;
	//	cout << "     Initial fragment" << G4endl;
	//	cout << theExcitedNucleus << G4endl << G4endl;

// 	G4LorentzVector MomentumTest = initialMomentum;
// 	G4cout << G4endl << G4endl << "Initial Momentum : "<< initialMomentum << G4endl;


	// DeExcite the nucleus 
	G4DynamicParticleVector * theDynamicPartVectorPtr = thePreCompoundModelPtr->DeExcite(theExcitedNucleus);
	
	G4double NofP = 0.0;
	G4double NofN = 0.0;
	
	for (G4int h = 0; h < theDynamicPartVectorPtr->entries(); h++) {
	  
// 	  G4cout << "Updating momentum ... " << MomentumTest << " - " << theDynamicPartVectorPtr->at(h)->Get4Momentum();
// 	  MomentumTest -= theDynamicPartVectorPtr->at(h)->Get4Momentum();
// 	  G4cout << " = " << MomentumTest <<theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() << G4endl;
	  
	  // 	theDynamicPartVectorPtr->at(h)->DumpInfo();
	  // 	G4cout << "PDGCharge: " << theDynamicPartVectorPtr->at(h)->GetDefinition()->GetPDGCharge() << G4endl;
	  

	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "proton") {
// 	    Protons->accumulate(theDynamicPartVectorPtr->at(h)->GetMass());
	    NofP++;
	    SumNP++;
	    ntupleP->column("P",theDynamicPartVectorPtr->at(h)->GetKineticEnergy());
	    ntupleP->column("angle",theDynamicPartVectorPtr->at(h)->GetMomentumDirection().angle(incidentProton));
	    ntupleP->dumpData();	    

	  }


	  
	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "neutron") {
// 	    Neutrons->accumulate(theDynamicPartVectorPtr->at(h)->GetMass());
	    SumNN++;
	    NofN++;
	    ntupleN->column("N",theDynamicPartVectorPtr->at(h)->GetKineticEnergy());
	    ntupleN->column("angle",theDynamicPartVectorPtr->at(h)->GetMomentumDirection().angle(incidentProton));
	    ntupleN->dumpData();	    
	  }
	  

	  
// 	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "alpha") {
// 	    NofP += 2;
// 	    SumNP += 2;
// 	    SumNN += 2;
// 	    NofN += 2;
// 	  }

// 	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "He3") {
// 	    NofP += 2;
// 	    SumNP += 2;
// 	    SumNN += 1;
// 	    NofN += 1;
// 	  }

// 	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "deuteron") {
// 	    NofP++;
// 	    SumNP++;
// 	    SumNN++;
// 	    NofN++;
// 	  }

// 	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "tritium") {
// 	    NofP++;
// 	    SumNP++;
// 	    SumNN += 2;
// 	    NofN += 2;
// 	  }



	  
// 	  if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() != "gamma") {
	    
// 	    ntuple->column("PDGEncoding",theDynamicPartVectorPtr->at(h)->GetDefinition()->GetPDGEncoding());
// 	    ntuple->column("DynMass",theDynamicPartVectorPtr->at(h)->GetMass());
// 	    ntuple->column("PDGMass",theDynamicPartVectorPtr->at(h)->GetDefinition()->GetPDGMass());
// 	    ntuple->column("Z",theDynamicPartVectorPtr->at(h)->GetDefinition()->GetPDGCharge());
// 	    ntuple->column("Px",theDynamicPartVectorPtr->at(h)->Get4Momentum().px());
// 	    ntuple->column("Py",theDynamicPartVectorPtr->at(h)->Get4Momentum().py());
// 	    ntuple->column("Pz",theDynamicPartVectorPtr->at(h)->Get4Momentum().pz());
// 	    ntuple->column("E",theDynamicPartVectorPtr->at(h)->Get4Momentum().e());
// 	    ntuple->dumpData();	    
// 	  }
	}

// 	ntuple2->column("P",NofP);
// 	ntuple2->column("N",NofN);
// 	ntuple2->dumpData();

// 	G4cout << "*************** PreEqTest *********************" << G4endl;
// 	G4cout << "Test Momentum: " << MomentumTest << G4endl;
// 	G4cout << "***********************************************" << G4endl;

      theDynamicPartVectorPtr->clearAndDestroy();
      delete theDynamicPartVectorPtr;
      }
      	
//       G4LorentzVector TestMomentum(initialMomentum);
//       for (G4int j=0; j < theFragVector->entries(); j++) {
// 	cout << theFragVector->at(j) << G4endl;

// 	// Test 4-momentum conservation 
// 	TestMomentum -= theFragVector->at(j)->GetMomentum();
//       }

//       cout << "******************" << G4endl;
      //       cout << "* Test Momentum = " << TestMomentum << G4endl;
//       cout << "******************" << G4endl;
      
      SumNP /= events;
      SumNN /= events;


      G4cout << "U = " << ExcitE[k] << "  Proton Multiplicity = " << SumNP << "  Neutron Multiplicity = " << SumNN << G4endl;



    }
delete [] ExcitE;
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

      G4double AtomicMass = G4ParticleTable::GetParticleTable()->GetIonTable()->GetIonMass(MyA,MyZ)/MeV;
      
      G4double MyExE = RandFlat::shoot(2.0*MyA,10.0*MyA);
      G4double MyPx = RandFlat::shoot(-2000,2000);
      G4double MyPy = RandFlat::shoot(-2000,2000);
      G4double MyPz = RandFlat::shoot(-2000,2000);

      G4int MyParticles = RandFlat::shoot(5);
      G4int MyHoles = RandFlat::shoot(5);
      G4int MyCharged = RandFlat::shoot(MyParticles);


      G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
      G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+
						(AtomicMass*MeV+MyExE*MeV)*
						(AtomicMass*MeV+MyExE*MeV))
				      );

    
      // put info about excited nucleus in fragment class
      G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);
      theExcitedNucleus.SetNumberOfExcitons(MyParticles+MyHoles);
      theExcitedNucleus.SetNumberOfHoles(MyHoles);
      theExcitedNucleus.SetNumberOfCharged(MyCharged);      


      cout << "Excited fragment: "<< G4endl;
      cout << theExcitedNucleus << G4endl;
      
      // DeExcite the nucleus 
      G4DynamicParticleVector * theDynamicPartVectorPtr = thePreCompoundModelPtr->DeExcite(theExcitedNucleus);      

      //      for (G4int h = 0; h < theDynamicPartVectorPtr->entries(); h++) theDynamicPartVectorPtr->at(h)->DumpInfo();


//       G4LorentzVector TestMomentum(initialMomentum);
//       for (G4int j=0; j < theFragVector->entries(); j++) {
// 	//	cout << theFragVector->at(j) << G4endl;

// 	// Test 4-momentum conservation 
// 	TestMomentum -= theFragVector->at(j)->GetMomentum();
//       }

//       cout << "******************" << G4endl;
//       cout << "* Test Momentum = " << TestMomentum << G4endl;
//       cout << "******************" << G4endl;



      theDynamicPartVectorPtr->clearAndDestroy();
      delete theDynamicPartVectorPtr;
      
    }
    
  }
  HbookManager->write();

  delete thePreCompoundModelPtr;
  delete theExcitationHandlerPtr;
  
    return 0;

}

