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
// gpaterno, October 2025
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "RunActionMessenger.hh"
#include "DetectorConstruction.hh"
#include "Run.hh"

#include "G4AnalysisManager.hh"
#include "G4RunManager.hh"
#include "G4Run.hh"
#include "G4RegionStore.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <iostream>
#include <fstream>
#include <string>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction():G4UserRunAction()
{
    G4RunManager::GetRunManager()->SetPrintProgress(1000);

    fMessenger = new RunActionMessenger(this);
    
    //Create the analysis manager  
    fAnalysisManager = G4AnalysisManager::Instance();      
    fAnalysisManager->SetDefaultFileType("root"); 
    fAnalysisManager->SetFileName(fFileName);
#ifdef G4MULTITHREADED
    fAnalysisManager->SetNtupleMerging(true);
#else
    fAnalysisManager->SetNtupleMerging(false);
#endif
    fAnalysisManager->SetVerboseLevel(0);
           
    //Creating the ntuple to score the particles impinging on the virtual screens
    fAnalysisManager->CreateNtuple("scoring_ntuple","virtual scoring screens");
    fAnalysisManager->CreateNtupleIColumn("screenID");
    fAnalysisManager->CreateNtupleSColumn("particle");
    fAnalysisManager->CreateNtupleDColumn("x");
    fAnalysisManager->CreateNtupleDColumn("y");
    fAnalysisManager->CreateNtupleDColumn("px");
    fAnalysisManager->CreateNtupleDColumn("py");
    fAnalysisManager->CreateNtupleDColumn("pz");
    fAnalysisManager->CreateNtupleDColumn("t");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();
        
    //Creating the ntuple to score Edep in the Radiatior
    fAnalysisManager->CreateNtuple("edep_rad","radiator");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();
    
    //Creating the ntuple to score Edep in the Converter
    fAnalysisManager->CreateNtuple("edep_conv","converter");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();
    
    //Creating the ntuple to score Edep in the spheres of the Granular Target
    fAnalysisManager->CreateNtuple("edep_spheres","granular target");
    fAnalysisManager->CreateNtupleIColumn("volumeID");
    fAnalysisManager->CreateNtupleDColumn("edep");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->FinishNtuple();

    //Creating the ntuple to score the particles exiting the crystals
    fAnalysisManager->CreateNtuple("scoring_ntuple2","particle leaving crystals");
    fAnalysisManager->CreateNtupleSColumn("particle");
    fAnalysisManager->CreateNtupleDColumn("x");
    fAnalysisManager->CreateNtupleDColumn("y");
    fAnalysisManager->CreateNtupleDColumn("z");
    fAnalysisManager->CreateNtupleDColumn("px");
    fAnalysisManager->CreateNtupleDColumn("py");
    fAnalysisManager->CreateNtupleDColumn("pz");
    fAnalysisManager->CreateNtupleDColumn("t");
    fAnalysisManager->CreateNtupleIColumn("eventID");
    fAnalysisManager->CreateNtupleIColumn("trackID");
    fAnalysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun() {return new Run;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
    //merging is only for root
    if (fFileName.find(".csv") != std::string::npos) {
        fAnalysisManager->SetNtupleMerging(false);
    }
       
    //open the output file         
    if (!fIsFileOpened) {
        fAnalysisManager->OpenFile(fFileName);
        fIsFileOpened = true;
    }

    //print begin message
    if (IsMaster()) {
        G4int NumberOfEventToBeProcessed = run->GetNumberOfEventToBeProcessed();
        
        G4cout
        << G4endl
        << "--------------------Begin of Global Run-----------------------" 
        << G4endl
        << "Number of events to be processed: " << NumberOfEventToBeProcessed       
        << G4endl
        << "--------------------------------------------------------------"
        << G4endl;
    }    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* run)
{
    //Write results on the root file
    if (fIsFileOpened) {
        fAnalysisManager->Write();
        fAnalysisManager->CloseFile();

        G4int threadID = G4Threading::G4GetThreadId();
        if (threadID > 0) {
            G4cout << "writing results of thread " << G4Threading::G4GetThreadId()
                   << " to root file: " << fFileName << G4endl;
        }
    }

    //Get the number of events
    G4int nofEvents = run->GetNumberOfEvent();
    if (nofEvents == 0) return;
                      
    //Get results from the Edep (in the Radiator Crystal) scorer
    const Run* aRun = static_cast<const Run*>(run);
    G4double Edep  = aRun->GetEdep();    
    G4double Edep2 = aRun->GetEdep2();     
    G4int nGoodEvents = aRun->GoodEvents();
    G4double stdEdep = std::sqrt(Edep2 - Edep*Edep/nGoodEvents); 
       
    //Print and write results to a text file
    if (IsMaster()) {    
        //detector construction instance    
        const DetectorConstruction* detectorConstruction = 
            static_cast<const DetectorConstruction*>(G4RunManager::GetRunManager()
                ->GetUserDetectorConstruction());

        G4String matName = 
            detectorConstruction->GetCrystalVolume()->GetMaterial()->GetName();            
        
        //set the significant digits to print
        G4cout.precision(8);
      
        //print the results on screen
        G4cout << G4endl
        << "--------------------End of Global Run-----------------------"
        << G4endl
        << " The run had " << nofEvents << " events";
        G4cout << G4endl
        << " Edep in " << "the Radiator Crystal" << " (made of "
        << matName << "): " << Edep/MeV
        << " +/- " << stdEdep/MeV << " MeV" << G4endl 
        << "------------------------------------------------------------" 
        << G4endl << G4endl;           
    }            
                   
    //end message
    G4cout
        << G4endl
        << " The run consisted of " << nofEvents << " particles" << G4endl
        << "------------------------------------------------------------"
        << G4endl << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetFileName(G4String filename) 
{
    if (filename != "") fFileName = filename;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

