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
// $Id: test15_hadronic_exerciser.cc,v 1.2 2009-10-20 07:39:16 kurasige Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

  G4int N = 10;
  if (argc > 1) {
    if (strcmp (argv[1], "large_N") == 0) {
      N = 200;
    }
  }

  G4UnitDefinition::BuildUnitsTable();

  G4cout <<
    "#"
    "\n# Auto-generated test input file for test15 hadronics."
    "\n#"
    "\n/control/verbose 2"
    "\n# /run/verbose 2"
    "\n/run/setCut 1 m"
    "\n/run/initialize"
    "\n/gun/direction 1 0 0";

  std::vector <G4String> particleNameList;
  particleNameList.push_back ("proton");

  std::vector <G4double> energyList;
  energyList.push_back (10 * GeV);

  std::vector <G4String> materialNameList;
  materialNameList.push_back ("Pb");
  materialNameList.push_back ("Al");
  materialNameList.push_back ("Air");
  materialNameList.push_back ("U");

  OutputCases (N, particleNameList, energyList, materialNameList);

  G4cout << G4endl;
}
