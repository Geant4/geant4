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
#include "globals.hh"
#include "G4ios.hh"
#include "G4qmdDummy.hh"

int main()
{

  G4String  InitialFile = "";

  G4cout << "Please enter filename of qMD data for initialization:" << G4endl;
  G4cin >> InitialFile;

	if (InitialFile == "") {
  	G4cout << "Initial file empty; using default file." << G4endl;
  }
	else {
  	G4cout << "Initial data file: " << InitialFile << G4endl;
	}

  G4qmdDummy theRun(InitialFile);
	
//  theRun.SetupFromFile();
//  theRun.justRun();
  G4cout << "   ...using inital file " << theRun.GetInputFile() << G4endl;

}



















