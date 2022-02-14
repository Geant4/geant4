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
// -------------------------------------------------------------------
// -------------------------------------------------------------------

#include "RunAction.hh"
#include "G4Run.hh"
#include "TrackingAction.hh"
#include "G4ParticleDefinition.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"
#include "G4Threading.hh"

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
{
 fFileName = "microelectronics";
 fpTrackingAction = 0;
 fInitialized = 0;
 fDebug = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* run)
{  
  // In this example, we considered that the same class was
  // used for both master and worker threads.
  // However, in case the run action is long,
  // for better code review, this practice is not recommanded.
  //
  
  if(isMaster) // WARNING : in sequential mode, isMaster == true    
    BeginMaster(run);    
  else 
    BeginWorker(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* run)
{
  if(isMaster)   
    EndMaster(run);
  else	
    EndWorker(run);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginMaster(const G4Run* run)
{
  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() == G4RunManager::sequentialRM);
  
  if(fDebug)
    {
      G4cout << "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << G4endl;
      if(!sequential)
	G4cout << "°°°°°°°°°°°°°°°° RunAction::BeginMaster" << G4endl;
      PrintRunInfo(run);
      G4cout << "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << G4endl;
    }
  
  if(sequential)
    {
      if(!fInitialized)	
	InitializeWorker(run);
      // Note: fpTrackingAction could be used as a flag for initialization instead      

      CreateHistogram();
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginWorker(const G4Run* run)
{
  if(fDebug)
    {
      G4cout << "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << G4endl;
      G4cout << "°°°°°°°°°°°°°°°° RunAction::BeginWorker" << G4endl;
      PrintRunInfo(run);
      G4cout << "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°" << G4endl;
    }
  if(!fInitialized)	
    InitializeWorker(run);
  
  CreateHistogram();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndMaster(const G4Run* run)
{
  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() 
		     == G4RunManager::sequentialRM);
  if(sequential)    
    EndWorker(run);    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndWorker(const G4Run* run)
{
  if(fDebug)
    {
      PrintRunInfo(run);
    }
  
  G4int nofEvents = run->GetNumberOfEvent();
  if ( nofEvents == 0 )
    {
      if(fDebug)
	{
	  G4cout << "°°°°°°°°°°°°°°°° NO EVENTS TREATED IN THIS RUN ==> LEAVING RunAction::EndOfRunAction "<< G4endl;
	}
      return;
    }
  
  ///////////////
  // Write Histo
  //
  WriteHistogram();

  ///////////////
  // Complete cleanup
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Clear();

  ///////////////
  // Printouts
  //
  std::map<const G4ParticleDefinition*, int>&
    particlesCreatedInWorld = fpTrackingAction->GetNParticlesCreatedInWorld();
  
  G4cout << "Number and type of particles created outside region \"Target\" :" << G4endl;
  
  PrintNParticles(particlesCreatedInWorld);
  
  G4cout << "_______________________" << G4endl;
  
  std::map<const G4ParticleDefinition*, int>&
    particlesCreatedInTarget = fpTrackingAction->GetNParticlesCreatedInTarget();
  
  G4cout << "Number and type of particles created in region \"Target\" :" << G4endl;
  
  PrintNParticles(particlesCreatedInTarget);
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::InitializeWorker(const G4Run*)
{
	if (fpTrackingAction == 0)
	{
		fpTrackingAction = (TrackingAction*) G4RunManager::GetRunManager()->GetUserTrackingAction();

		if(fpTrackingAction == 0 && isMaster == false)
		{
			G4ExceptionDescription exDescrption ;
			exDescrption << "fpTrackingAction is a null pointer. Has it been correctly initialized ?";
			G4Exception("RunAction::BeginOfRunAction","RunAction001",FatalException, exDescrption);
		}
	}

	fInitialized = true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateHistogram()
{
	// Book histograms, ntuple

	// Create analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in Analysis.hh

	G4cout << "##### Create analysis manager " << "  " << this << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");

	G4cout << "Using " << analysisManager->GetType() << " analysis manager" << G4endl;

	// Create directories

	//analysisManager->SetHistoDirectoryName("histograms");
	//analysisManager->SetNtupleDirectoryName("ntuple");
	analysisManager->SetVerboseLevel(1);

	// Open an output file
	
	analysisManager->OpenFile(fFileName);

	// Creating ntuple

	analysisManager->CreateNtuple("microelectronics", "physics");
	analysisManager->CreateNtupleDColumn("flagParticle");
	analysisManager->CreateNtupleDColumn("flagProcess");
	analysisManager->CreateNtupleDColumn("x");
	analysisManager->CreateNtupleDColumn("y");
	analysisManager->CreateNtupleDColumn("z");
	analysisManager->CreateNtupleDColumn("totalEnergyDeposit");
	analysisManager->CreateNtupleDColumn("stepLength");
	analysisManager->CreateNtupleDColumn("kineticEnergyDifference");
	analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::WriteHistogram()
{
	// print histogram statistics
	//
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

	// save histograms
	//
	analysisManager->Write();
	analysisManager->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::PrintRunInfo(const G4Run* run)
{
	G4cout << "°°°°°°°°°°°°°°°° Run is = " << run->GetRunID() << G4endl;
	G4cout << "°°°°°°°°°°°°°°°° Run type is = " << G4RunManager::GetRunManager()->GetRunManagerType() << G4endl;
	G4cout << "°°°°°°°°°°°°°°°° Event processed = " << run->GetNumberOfEventToBeProcessed() << G4endl;
	G4cout << "°°°°°°°°°°°°°°°° N° Event = " << run->GetNumberOfEvent() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container)
{
    std::map<const G4ParticleDefinition*, int>::iterator it;
    for(it = container.begin() ;
        it != container.end(); it ++)
    {
        G4cout << "N " << it->first->GetParticleName() << " : " << it->second << G4endl;
    }
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
