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
/// \file exoticphysics/monopole/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id: RunAction.cc 107534 2017-11-21 13:13:36Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "RunActionMessenger.hh" 
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"

#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "Randomize.hh"

#include "G4EmCalculator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
  :fDetector(det),fKinematic(kin)
{
  fMessenger = new RunActionMessenger(this);
  fBinLength = 5 * CLHEP::mm; 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  analysisManager->SetFileName("monopole");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
#ifdef G4MULTITHREADED
  if(isMaster) delete fKinematic;
#endif
  delete fMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector,fKinematic);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //histograms
  //        
  Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
  // print Run summary
  //
  if (isMaster) fRun->EndOfRun(fBinLength);    
      
  // save histograms
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {    
    analysisManager->Write();
    analysisManager->CloseFile();
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetBinSize(G4double size)
{ 
  fBinLength = size;
  if(fBinLength > fDetector->GetMaxStepSize()) { 
    fBinLength = fDetector->GetMaxStepSize();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Book()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  analysisManager->SetFirstHistoId(1);   

  G4double length = fDetector->GetAbsorSizeX();
  G4int nbBins = G4lrint(length / fBinLength);

  // Create histograms
  analysisManager->CreateH1("h1","Edep (MeV/mm) along absorber (mm)", 
                            nbBins, 0, length);
  analysisManager->CreateH1("h2","DEDX (MeV/mm) of proton", 100, -3., 7.);
  analysisManager->CreateH1("h3","DEDX (MeV/mm) of monopole", 100, -3., 7.);
  analysisManager->CreateH1("h4","Range(mm) of proton", 100, -3., 7., "mm");
  analysisManager->CreateH1("h5","Range(mm) of monopole", 100, -3., 7., "mm");
  analysisManager->OpenFile(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
