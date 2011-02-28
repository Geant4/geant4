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
//

#include "globals.hh"
#include "G4UnitsTable.hh"
#include <vector>

static void OutputCases
(G4int N,  // Number of events per case.
 const std::vector <G4String> & particleNameList,
 const std::vector <G4double> & energyList,
 const std::vector <G4String> & materialNameList) {

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

  G4int N = 100;
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
    "\n/run/setCut 1 m"
    "\n/run/initialize"
    "\n/gun/direction 1 0 0";

  std::vector <G4String> particleNameList;
  particleNameList.push_back ("proton");
  particleNameList.push_back ("neutron");
  particleNameList.push_back ("pi+");
  particleNameList.push_back ("pi-");
  particleNameList.push_back ("kaon+");
  particleNameList.push_back ("kaon-");
  particleNameList.push_back ("kaon0S");
  particleNameList.push_back ("kaon0L");

  std::vector <G4double> energyList;
  energyList.push_back (100 * GeV);

  std::vector <G4String> materialNameList;
  materialNameList.push_back ("Pb");
  materialNameList.push_back ("Al");
  materialNameList.push_back ("Air");

  OutputCases (N, particleNameList, energyList, materialNameList);

// anti - particles testing... 
  std::vector <G4String> anti_particleNameList;
  anti_particleNameList.push_back ("anti_proton");
  anti_particleNameList.push_back ("anti_neutron");
  anti_particleNameList.push_back ("anti_deuteron");
  anti_particleNameList.push_back ("anti_triton");
  anti_particleNameList.push_back ("anti_He3");
  anti_particleNameList.push_back ("anti_alpha");

  std::vector <G4double> anti_energyList;
  anti_energyList.push_back (10 * GeV);

  OutputCases (5, anti_particleNameList, anti_energyList, materialNameList);


  G4cout << G4endl;
}
