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


#include "G4Evaporation.hh"
#include "G4NucleiProperties.hh"
#include "G4Fragment.hh"
#include "G4ThreeVector.hh"

int main()
{
  cout << G4endl << "Evaporation Test program" << G4endl;


  // HBOOK initialization
  HepTupleManager* HbookManager;
  HbookManager = new HBookFile("evaptest.hbook", 58);
  // Book a ntuple
  HepTuple* ntupleP = HbookManager->ntuple("Protons");;
  HepTuple* ntupleN =  HbookManager->ntuple("Neutrons");;
 
  G4Evaporation * theEvaporationPtr = new G4Evaporation;

  // Excited nucleus
  G4Fragment theExcitedNucleus;

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
    G4double AtomicMass = G4NucleiProperties::GetAtomicMass(MyA,MyZ);
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
    G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+
					      MeV*MeV*(AtomicMass+MyExE)*(AtomicMass+MyExE)));
    
    // put info about excited nucleus in fragment class
    theExcitedNucleus.SetA(MyA);
    theExcitedNucleus.SetZ(MyZ);
    theExcitedNucleus.SetExcitationEnergy(MyExE*MeV); // in G4 units
    theExcitedNucleus.SetMomentum(initialMomentum);
    // end excited nucleus

    for (G4int j = 0; j < events; j++) {
      // DeExcite the nucleus (it returns a vector of dynamical particles)

      G4FragmentVector * theFragVector = theEvaporationPtr->BreakItUp(theExcitedNucleus);


      for (G4int h = 0; h < theFragVector->entries(); h++) {
	if (theFragVector->at(h)->GetA() == 1 && theFragVector->at(h)->GetZ() == 0) {
	  ntupleN->column("N",theFragVector->at(h)->GetMomentum().e()-G4NucleiProperties::GetAtomicMass(1,0));
	  ntupleN->dumpData();	    
	  
	} else if (theFragVector->at(h)->GetA() == 1 && theFragVector->at(h)->GetZ() == 1) {
	  ntupleP->column("P",theFragVector->at(h)->GetMomentum().e()-G4NucleiProperties::GetAtomicMass(1,1));
	  ntupleP->dumpData();	    
	}
      }




      theFragVector->clearAndDestroy();
      delete theFragVector;
    }
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
      
      G4double MyExE = RandFlat::shoot(2000);
      G4double MyPx = RandFlat::shoot(-2000,2000);
      G4double MyPy = RandFlat::shoot(-2000,2000);
      G4double MyPz = RandFlat::shoot(-2000,2000);

      G4ThreeVector triV(MyPx*MeV,MyPy*MeV,MyPz*MeV);
      G4LorentzVector initialMomentum(triV,sqrt(triV.mag2()+AtomicMass*AtomicMass));
    
      // put info about excited nucleus in fragment class
      theExcitedNucleus.SetA(MyA);
      theExcitedNucleus.SetZ(MyZ);
      theExcitedNucleus.SetExcitationEnergy(MyExE*MeV); 
      theExcitedNucleus.SetMomentum(initialMomentum);
      // end excited nucleus
      
      
      cout << "Excited fragment: "<< G4endl;
      cout << theExcitedNucleus << G4endl;
      
      G4FragmentVector * theFragVector = theEvaporationPtr->BreakItUp(theExcitedNucleus);
      
      theFragVector->clearAndDestroy();
      delete theFragVector;
      
    }
    
  }
  HbookManager->write();

  delete theEvaporationPtr;
  
  return 0;

}
