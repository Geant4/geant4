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
#include "G4ios.hh"
#include "globals.hh"

#include "CLHEP/Hist/TupleManager.h"
#include "CLHEP/Hist/HBookFile.h"
#include "CLHEP/Hist/Histogram.h"
#include "CLHEP/Hist/Tuple.h"

#include "G4ExcitationHandler.hh"
#include "G4NucleiProperties.hh"
#include "G4DynamicParticle.hh"
#include "G4DynamicParticleVector.hh"
#include "G4Fragment.hh"
#include "G4ThreeVector.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4Electron.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"


// #include "CLHEP/String/Strings.h"


// void prueba(void) 
// {
//   HepString fileName;
  
//   for (G4int i = 0; i < 1000; i++) {
//     HepString nameZ = "/z" + HepString(i);
//     HepString nameA = ".a" + HepString(i);
//     fileName = nameZ + nameA;
//   }

//   return;
// }


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

  G4cout << G4endl << "De-Excitation Test program" << G4endl;

  // HBOOK initialization
  HepTupleManager* HbookManager;
  HbookManager = new HBookFile("exhandler.hbook", 58);
  // Book a ntuple
  HepTuple* ntupleP = HbookManager->ntuple("Protons");
  HepTuple* ntupleN =  HbookManager->ntuple("Neutrons");



  G4ParticleTable * theParticleTable =   G4ParticleTable::GetParticleTable();

  
  G4ExcitationHandler * theHandler = new G4ExcitationHandler;


  cout << G4endl << G4endl << G4endl;

   
  cout << "A = ";
  G4int MyA;
  G4cin >> MyA;
  cout << "Z = ";
  G4int MyZ;
  G4cin >> MyZ;
  cout << "ExcitationEnergy (MeV) = ";
  G4double MyExE;
  G4cin >> MyExE;
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
  

  G4int events = 1;
  cout << "Iterations: ";
  G4cin >> events;
  
    
  G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
  G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
    
  // put info about excited nucleus in fragment class

  G4Fragment theExcitedNucleus(MyA,MyZ,initialMomentum);


  for (G4int j = 0; j < events; j++) {

    G4cout << "Iteration: " << j+1 << G4endl;
    // DeExcite the nucleus (it returns a vector of dynamical particles)
    G4DynamicParticleVector * theDynamicPartVectorPtr  = theHandler->BreakItUp(theExcitedNucleus);
    
    for (G4int h = 0; h < theDynamicPartVectorPtr->entries(); h++) {

      if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "proton") {
	ntupleP->column("P",theDynamicPartVectorPtr->at(h)->GetKineticEnergy());
	ntupleP->dumpData();	      
      }    
    
      if (theDynamicPartVectorPtr->at(h)->GetDefinition()->GetParticleName() == "neutron") {
	ntupleN->column("N",theDynamicPartVectorPtr->at(h)->GetKineticEnergy());
	ntupleN->dumpData();	    
      }
    }
      
    theDynamicPartVectorPtr->clearAndDestroy();
    delete theDynamicPartVectorPtr;
  }
  
  HbookManager->write();
  
  delete theHandler;
  
  //  prueba();
  
  return EXIT_SUCCESS;

}

