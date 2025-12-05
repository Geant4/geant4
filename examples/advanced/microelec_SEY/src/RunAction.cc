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
#include "PrimaryGeneratorAction.hh"
#include "MicroElecRun.hh"
#include "G4UnitsTable.hh"
#include "g4csv_defs.hh"
#include "G4SystemOfUnits.hh"

void PrintNParticles(std::map<const G4ParticleDefinition*, int>& container);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::RunAction()
{
 fFileName = "microelectronics";
 SEYfileName = "SEY";
 SEYntupleName = "SEYntuple";
 fpTrackingAction = 0;
 fInitialized = 0;
 NbOfIncidentPart = 0;
 EnergyOfIncidentPart = 999;
 fDebug = false;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4Run* RunAction::GenerateRun()
{
	fRun = new MicroElecRun();
	return fRun;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::BeginOfRunAction(const G4Run* run)
{  
  // In this example, we considered that the same class was
  // used for both master and worker threads.
  // However, in case the run action is long,
  // for better code review, this practice is not recommanded.
  //
	if (isMaster) {// WARNING : in sequential mode, isMaster == true    
		//PrimaryGeneratorAction* PGA = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
		//EnergyOfIncidentPart = PGA->GetEnergy();
		//EnergyOfIncidentPart = PGA->GetParticleGPS()->GetParticleEnergy();
		BeginMaster(run);
	}
  else 
    BeginWorker(run);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndOfRunAction(const G4Run* run)
{

	if (isMaster) {
		EndMaster(run);
		G4cout << "End of run action : " << EnergyOfIncidentPart << G4endl;
	}
	else {
	PrimaryGeneratorAction* PGA = (PrimaryGeneratorAction*)G4RunManager::GetRunManager()->GetUserPrimaryGeneratorAction();
	EnergyOfIncidentPart = PGA->GetParticleGPS()->GetCurrentSource()->GetEneDist()->GetMonoEnergy();
	EndWorker(run);
	}

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

      //CreateHistogram();
    }
  // modif CI 1/3/2022
  else {
	  if (run->GetRunID() == 0) { CreateSEYHistogram(); }
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
  // modif CI 4/3/2022
  //CreateHistogram();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::EndMaster(const G4Run* run)
{
  bool sequential = (G4RunManager::GetRunManager()->GetRunManagerType() 
		     == G4RunManager::sequentialRM);
  if (sequential) {
	  EndWorker(run);
  }
  else {
	  
	  
	  
	  if (run->GetRunID() == 0) {
		  NbOfIncidentPart  = run->GetNumberOfEvent();
	  }
	  

	  MicroElecRun* localRun = (MicroElecRun*)(run);

	  EnergyOfIncidentPart = localRun->GetElecEneIncPart();
	  G4double compteurTot = localRun->GetElecTotaScorer();
	  G4double compteurPrimaire = localRun->GetElecPrimScorer();
	  G4double compteurSec = localRun->GetElecSecoScorer();
	  G4double compteurSup50 = localRun->GetElecSup50Scorer();

	  G4cout << "ooo End Master ooo, Run : " << run->GetRunID() << " | Energy = " << EnergyOfIncidentPart / eV <<
		  " eV | Number of inci Part. = " << NbOfIncidentPart << " | Tot = " << compteurTot <<
		  " | Primaires = " << compteurPrimaire << " | Sec = " << compteurSec << G4endl <<
		  " | Taux d'emission = " << compteurTot/ NbOfIncidentPart << G4endl;

	  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	  analysisManager->FillNtupleDColumn(0, EnergyOfIncidentPart/eV);
	  analysisManager->FillNtupleDColumn(1, compteurTot/ NbOfIncidentPart);
	  analysisManager->FillNtupleDColumn(2, compteurSec/ NbOfIncidentPart);
	  analysisManager->FillNtupleDColumn(3, compteurPrimaire/ NbOfIncidentPart);
	  analysisManager->FillNtupleDColumn(4, compteurSup50/ NbOfIncidentPart);
	  analysisManager->AddNtupleRow();//*/
  }

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
  // Complete cleanup
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Clear();

  ///////////////
  // Printouts
  //
  std::map<const G4ParticleDefinition*, int>&
    particlesCreatedInWorld = fpTrackingAction->GetNParticlesCreatedInWorld();
  std::map<const G4ParticleDefinition*, int>&
	  particlesCreatedInTarget = fpTrackingAction->GetNParticlesCreatedInTarget();

  G4cout << "__________ooo End Worker ooo begin_____________" << G4endl;
  G4cout << " Number and type of particles created outside region \"Target\" :";
  PrintNParticles(particlesCreatedInWorld);
  G4cout << "ooo End Worker ooo, Number and type of particles created in region \"Target\" :";
  PrintNParticles(particlesCreatedInTarget);
  G4cout << "__________ooo End Worker ooo end_____________" << G4endl;
 

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


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::SetSEYFileName(G4String& name)
{
	SEYfileName = name;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::CreateSEYHistogram()
{
	// Book histograms, ntuple

	// Create analysis manager
	// The choice of analysis technology is done via selection of a namespace
	// in Analysis.hh

	G4cout << "##### Create analysis manager " << "  " << this << G4endl;
	G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
	analysisManager->SetDefaultFileType("csv");

	G4cout << "Using " << analysisManager->GetType() << " analysis manager" << G4endl;

	// Create directories

	//analysisManager->SetHistoDirectoryName("histograms");
	//analysisManager->SetNtupleDirectoryName("ntuple");
	analysisManager->SetVerboseLevel(1);

	// Open an output file

	analysisManager->OpenFile(SEYfileName);
	// Creating ntuple
	//									File name ext.   tupple info
	analysisManager->CreateNtuple("data", "Sec. Ele. Emission");
	analysisManager->CreateNtupleDColumn("Initial energy (eV)");
	analysisManager->CreateNtupleDColumn("TEEY");
	analysisManager->CreateNtupleDColumn("SEEY");
	analysisManager->CreateNtupleDColumn("BEEY");
	analysisManager->CreateNtupleDColumn("EEY>50eV");
	analysisManager->FinishNtuple();
	
	
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void RunAction::WriteSEYHistogram()
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
