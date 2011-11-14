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
// $Id: RunAction.cc,v 1.1 2009/03/21 18:37:27 vnivanch Exp $
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#ifdef G4_USE_ROOT
#include "ROOTAnalysis.hh"
#endif

#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Version.hh"

#include "PhysicsList.hh"
#include "Randomize.hh"

#include <time.h>    //needed on Windows for time()

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(PhysicsList* pl) : 
  outFile(0),fPL(pl)
{
  primaryEnergy = 1.0*keV;
  runID = -1;
  fRandomSeed = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  if (outFile)
    {
      outFile->close();
      delete outFile;
    }

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  runID = aRun->GetRunID();

  //For Run=0, check the random seed
  if (!runID)
    {
      if (fRandomSeed < 0) // not set by anyone!
	{
	  fRandomSeed = time(0);
	  G4cout << "Choosing random seed according to time(0) " << G4endl;	  
	}
      else
	G4cout << "Constant random seed: " << fRandomSeed << G4endl;

      CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
      CLHEP::HepRandom::setTheSeed(fRandomSeed);  
      CLHEP::HepRandom::showEngineStatus();      
    }

  if (!outFile)
    {
      outFile = new std::ofstream();
      //Use the version name for the dat file
      G4String filename;
      //remove blank spaces
      for (size_t i = 0; i < G4Version.length(); i++)
        if (G4Version[i] != ' ' && G4Version[i] != '$') filename += G4Version[i];
      //check if there is "Name" in front, if so remove 5 chars      
      if (filename.substr(0,4) == "Name")
	{
	  G4String vs = filename;
	  filename = vs.substr(5,vs.size()-5);
	}
      filename += "_";
      filename += fPL->GetEmName();
      filename += ".dat";
      outFile->open(filename.c_str());
    }

  PrimaryGeneratorAction* primGen = (PrimaryGeneratorAction*) 
    G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();

  primaryEnergy = primGen->GetPrimaryEnergy();
  
  G4cout << "Primary energy: " << primaryEnergy/keV << " keV " << G4endl;

#ifdef G4_USE_ROOT
  ROOTAnalysis::getInstance()->BookNewHistogram(aRun->GetRunID(),primaryEnergy);
#endif
 
  counter = 0;
  counterTot = 0;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4double nEvents = (G4double) (aRun->GetNumberOfEvent());
  if (nEvents)
    {           
      if (outFile->is_open())
	{
	  G4double dCounter = (G4double) counter;
	  (*outFile) << primaryEnergy/keV << " " << counter/nEvents << " " << 
	    std::sqrt(dCounter)/nEvents << G4endl;
	}
      else
	{
	  G4Exception("RunAction::EndOfRunAction()","tst67_05",
		      FatalException,"Unable to open output file");
	  return;
	}

#ifdef G4_USE_ROOT
      ROOTAnalysis::getInstance()->EndOfRun(runID,primaryEnergy,
					    (G4int) nEvents,
					    counter,counterTot);
#endif
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EventEnergy(G4double ene)
{
  if (ene < 0.1*eV) 
    return;
  counterTot++;

#ifdef G4_USE_ROOT
  ROOTAnalysis::getInstance()->AddEventEnergy(runID,ene);
#endif
  if (std::fabs(primaryEnergy - ene) < 0.001*primaryEnergy) //within 0.1%
    counter++;
}


