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
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
//

#include "globals.hh"
#include "G4UnitsTable.hh"
#include "g4std/vector"

static void OutputCases
(G4int N,  // Number of events per case.
 const G4std::vector <G4String> & particleNameList,
 const G4std::vector <G4double> & energyList,
 const G4std::vector <G4String> & materialNameList) {

  for (size_t iMaterial = 0;
       iMaterial < materialNameList.size ();
       iMaterial++) {
    for (size_t iEnergy = 0;
	 iEnergy < energyList.size ();
	 iEnergy++) {
      for (size_t iParticle = 0;
	   iParticle < particleNameList.size ();
	   iParticle++) {

	G4cout
	  << "\n#\n# " << particleNameList [iParticle]
	  << " at " << G4BestUnit (energyList [iEnergy], "Energy")
	  << " in " << materialNameList [iMaterial]
	  << "\n#";

	G4cout
	  << "\n/gun/particle " << particleNameList [iParticle]
	  << "\n/gun/energy " <<  G4BestUnit (energyList [iEnergy], "Energy")
	  << "\n/mydet/SelectMaterial " << materialNameList [iMaterial]
	  << "\n/run/beamOn " << N;

      }
    }
  }
}


int main (int argc, char** argv) {

  G4int N = 200;
  if (argc > 1) {
    if (strcmp (argv[1], "large_N") == 0) {
      N = 2000;
    }
  }

  G4UnitDefinition::BuildUnitsTable();

  G4cout <<
    "#"
    "\n# Auto-generated test input file for test12 hadronics."
    "\n#"
    "\n/control/verbose 2"
    "\n# /run/verbose 2"
    "\n/run/particle/setCut 1 m"
    "\n/run/initialize"
    "\n/gun/direction 1 0 0";

  G4std::vector <G4String> particleNameList;
  particleNameList.push_back ("proton");
  particleNameList.push_back ("neutron");
  particleNameList.push_back ("pi+");
  particleNameList.push_back ("pi-");
  particleNameList.push_back ("kaon+");
  particleNameList.push_back ("kaon-");
  particleNameList.push_back ("kaon0S");
  particleNameList.push_back ("kaon0L");

  G4std::vector <G4double> energyList;
  energyList.push_back (100 * GeV);

  G4std::vector <G4String> materialNameList;
  materialNameList.push_back ("Pb");
  materialNameList.push_back ("Al");
  materialNameList.push_back ("Air");

  OutputCases (N, particleNameList, energyList, materialNameList);

  G4cout << G4endl;
}
