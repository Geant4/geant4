// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: test02.hadronic.exerciser.cc,v 1.5 1999-12-15 14:54:36 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "g4rw/tvordvec.h"

static void OutputCases
(const G4RWTValOrderedVector <G4String> & particleNameList,
 const G4RWTValOrderedVector <G4double> & energyList,
 const G4RWTValOrderedVector <G4String> & materialNameList) {

  for (int iMaterial = 0;
       iMaterial < materialNameList.entries ();
       iMaterial++) {
    for (int iEnergy = 0;
	 iEnergy < energyList.entries ();
	 iEnergy++) {
      for (int iParticle = 0;
	   iParticle < particleNameList.entries ();
	   iParticle++) {

	cout << "\n#\n# " << particleNameList [iParticle]
	     << " at " << G4BestUnit (energyList [iEnergy], "Energy")
	     << " in " << materialNameList [iMaterial]
	     << "\n#";

	cout
	  << "\n/gun/particle " << particleNameList [iParticle]
	  << "\n/gun/energy " <<  G4BestUnit (energyList [iEnergy], "Energy")
	  << "\n/mydet/SelectMaterial " << materialNameList [iMaterial]
	  << "\n/run/beamOn 1";

      }
    }
  }
}


int main (int argc, char** argv) {

  G4UnitDefinition::BuildUnitsTable();

  cout <<
    "#"
    "\n# Auto-generated test input file for test02 hadronics."
    "\n#"
    "\n/control/verbose 2"
    "\n# /run/verbose 2"
    "\n/run/initialize"
    "\n/gun/direction 1 0 0";

  G4RWTValOrderedVector <G4String> particleNameList;
  particleNameList.append ("proton");
  particleNameList.append ("pi+");
  particleNameList.append ("pi-");
  particleNameList.append ("kaon+");
  particleNameList.append ("kaon-");
  particleNameList.append ("kaon0S");
  particleNameList.append ("kaon0L");

  G4RWTValOrderedVector <G4double> energyList;
  energyList.append (100 * GeV);

  G4RWTValOrderedVector <G4String> materialNameList;
  materialNameList.append ("Pb");
  materialNameList.append ("Al");
  materialNameList.append ("Air");

  OutputCases (particleNameList, energyList, materialNameList);

  particleNameList.clear ();
  particleNameList.append ("proton");
  particleNameList.append ("neutron");

  energyList.clear ();
  energyList.append (500 * MeV);

  materialNameList.clear ();
  materialNameList.append ("Pb");
  materialNameList.append ("Al");
  materialNameList.append ("Air");

  OutputCases (particleNameList, energyList, materialNameList);

  particleNameList.clear ();
  particleNameList.append ("neutron");
  particleNameList.append ("triton");
  particleNameList.append ("deuteron");
  particleNameList.append ("alpha");

  energyList.clear ();
  energyList.append (100 * MeV);

  materialNameList.clear ();
  materialNameList.append ("Pb");
  materialNameList.append ("Al");
  materialNameList.append ("Air");

  OutputCases (particleNameList, energyList, materialNameList);

  cout << G4endl;
}
