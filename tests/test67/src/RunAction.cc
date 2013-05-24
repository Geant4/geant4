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
// $Id$
// GEANT4 tag $Name:  $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include <time.h>    //needed on Windows for time()

#include "RunAction.hh"

#ifdef G4_USE_ROOT
#include "ROOTAnalysis.hh"
#endif

#include "G4SystemOfUnits.hh"

#include "PrimaryGeneratorAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4Version.hh"

#include "PhysicsList.hh"
#include "Randomize.hh"
#include "Run.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction() : 
  outFile(0)
{
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
G4Run* RunAction::GenerateRun()
{
  //this is needed for master and slaves.
  return new Run();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  const Run* theRun = static_cast<const Run*> (aRun);
  G4int runID = theRun->GetRunID();

#ifdef G4_USE_ROOT 
  //If it is a slave OR the sequential, book histo
  if (!IsMaster() || 
      G4RunManager::GetRunManager()->GetRunManagerType()==G4RunManager::sequentialRM)
    {      
      //create histos, if there are not already there for this run
      if (!(ROOTAnalysis::getInstance()->AreHistoCreated()))
	ROOTAnalysis::getInstance()->BookNewHistogram(theRun->GetRunID(),
						      theRun->GetPrimaryEnergy());
    }
#endif

  //Master or sequential
  if (IsMaster()) 
    G4cout << "ooo Run " << theRun->GetRunID() << " starts (global)." << G4endl;
  else //it is a slave, do nothing else
    {
      G4cout << "ooo Run " << theRun->GetRunID() << " starts on slave." << G4endl;
      return;
    }
    
  //For Run=0, check the random seed in the master
  if (!runID)
    {
      if (fRandomSeed < 0) // not set by anyone!
	{
	  fRandomSeed = time(0);
	  G4cout << "Choosing random seed according to time(0) " << G4endl;	  
	}
      else
	G4cout << "Constant random seed: " << fRandomSeed << G4endl;

      G4Random::setTheEngine(new CLHEP::RanecuEngine);
      G4Random::setTheSeed(fRandomSeed);  
      G4Random::showEngineStatus();      
    }

  if (!outFile)
    {
      //Master RunAction needs the physics list info
      const PhysicsList* physicsList =
	dynamic_cast<const PhysicsList*>
	(G4RunManager::GetRunManager()->GetUserPhysicsList());
      
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
      filename += physicsList->GetEmName();
      filename += ".dat";
      outFile->open(filename.c_str());
    }

#ifdef G4_USE_ROOT
  ROOTAnalysis::getInstance()->ResetHistoForNewRun(); //histo created by the first worker
#endif
   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  const Run* theRun = static_cast<const Run*> (aRun);
  
  //Master or sequential
  if (IsMaster() || 
      G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::sequentialRM)
    G4cout << "Global results with " << theRun->GetNumberOfEvent() << 
      " events : " << theRun->GetCounter() << G4endl;
  else //it is a slave, do nothing
    {
      G4cout << "Local results with " << theRun->GetNumberOfEvent() << 
	" events : " << theRun->GetCounter() << G4endl;
      return;
    }

  G4double nEvents = (G4double) (theRun->GetNumberOfEvent());

  if (nEvents)
    {           
      if (outFile->is_open())
	{
	  G4double dCounter = (G4double) theRun->GetCounter();
	  (*outFile) << theRun->GetPrimaryEnergy()/keV << " " << theRun->GetCounter()/nEvents 
		     << " " << 
	    std::sqrt(dCounter)/nEvents << G4endl;
	}
      else
	{
	  G4Exception("RunAction::EndOfRunAction()","tst67_05",
		      FatalException,"Unable to open output file");
	  return;
	}

#ifdef G4_USE_ROOT
      ROOTAnalysis::getInstance()->EndOfRun(theRun->GetRunID(),
					    theRun->GetPrimaryEnergy(),
					    (G4int) nEvents,
					    theRun->GetCounter(),
					    theRun->GetCounterTot());
#endif
    }
}
