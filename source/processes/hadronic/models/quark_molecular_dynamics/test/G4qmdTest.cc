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
#include "G4qmd.hh"
#include "G4KineticTrack.hh"
#include "G4KineticTrackVector.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

int main()
{
	G4KineticTrackVector * ResultingHadrons = new G4KineticTrackVector();

  G4String & InitialFile = "";
  G4double OutputStepping = 0.5;

  G4cout << "Please enter filename of qMD data for initialization:" << G4endl;
  G4cin >> InitialFile;

  G4cout << "-----------------------------------------" << G4endl;
  G4cout << "Please enter time steps of output:" << G4endl;
  G4cin >> OutputStepping;

	if (InitialFile == "") {
  	G4cout << "Initial file empty; using default file." << G4endl;
  }
	else {
  	G4cout << "Initial data file: " << InitialFile << G4endl;
	}



  G4qmd theRun(InitialFile);

  theRun.SetupFromFile();
  theRun.SetOutputTimestep(OutputStepping);

  ResultingHadrons = theRun.TheHadrons();


}



















