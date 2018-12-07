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
/// \file electromagnetic/TestEm2/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "Run.hh"
#include "EmAcceptance.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
 :G4UserRunAction(),
  fDet(det),fKin(kin),fAnalysisManager(nullptr),fRun(nullptr),
  fVerbose(0), fEdeptrue(1.), fRmstrue(1.), fLimittrue(DBL_MAX)
{
  fRunMessenger = new RunActionMessenger(this);
  fHistoName[0] = "testem2";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BookHisto()
{
  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  //

  fAnalysisManager = G4AnalysisManager::Instance();
    
  // Open an output file
  //
  ///fHistoName[0] = "testem2";

  fAnalysisManager->OpenFile(fHistoName[0]); 

  fAnalysisManager->SetVerboseLevel(1);
  G4String extension = fAnalysisManager->GetFileType();
  fHistoName[1] = fHistoName[0] + "." + extension;  

  // Creating histograms
  //
  G4double Ekin = fKin->GetParticleGun()->GetParticleEnergy();
  G4int    nLbin = fDet->GetnLtot();
  G4int    nRbin = fDet->GetnRtot();
  G4double dLradl = fDet->GetdLradl();
  G4double dRradl = fDet->GetdRradl();
  
  fAnalysisManager->SetFirstHistoId(1);   
  fAnalysisManager->CreateH1( "h1","total energy deposit(percent of Einc)",
                                  110,0.,110.);

  fAnalysisManager->CreateH1( "h2","total charged tracklength (radl)",
                                  110,0.,110.*Ekin/GeV);

  fAnalysisManager->CreateH1( "h3","total neutral tracklength (radl)",
                                  110,0.,1100.*Ekin/GeV);

  fAnalysisManager->CreateH1( "h4","longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);
                                    
  fAnalysisManager->CreateP1( "p4","longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl, 0., 1000.);
                                    
  fAnalysisManager->CreateH1( "h5","rms on longit Edep (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);

  G4double Zmin=0.5*dLradl, Zmax=Zmin+nLbin*dLradl;
  fAnalysisManager->CreateH1( "h6","cumul longit energy dep (% of E inc)",
                                  nLbin,Zmin,Zmax);                          
                                    
  fAnalysisManager->CreateH1( "h7","rms on cumul longit Edep (% of E inc)",
                                  nLbin,Zmin,Zmax);

  fAnalysisManager->CreateH1( "h8","radial energy profile (% of E inc)",
                                  nRbin,0.,nRbin*dRradl);
                                  
  fAnalysisManager->CreateP1( "p8","radial energy profile (% of E inc)",
                                  nRbin,0.,nRbin*dRradl, 0., 1000.);
                                  
  fAnalysisManager->CreateH1( "h9","rms on radial Edep (% of E inc)",
                                  nRbin,0.,nRbin*dRradl);            

  G4double Rmin=0.5*dRradl, Rmax=Rmin+nRbin*dRradl;
  fAnalysisManager->CreateH1("h10","cumul radial energy dep (% of E inc)",
                                  nRbin,Rmin,Rmax);

  fAnalysisManager->CreateH1("h11","rms on cumul radial Edep (% of E inc)",
                                  nRbin,Rmin,Rmax);                    
                                    
  G4cout << "\n----> Histogram file is opened in " << fHistoName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetHistoName(G4String& val)
{
  fHistoName[0] = val;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDet, fKin); 
  fRun->SetVerbose(fVerbose);
  return fRun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  // show Rndm status
  if (isMaster) G4Random::showEngineStatus(); 

  //histograms
  //
  BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
 //compute and print statistic
 //
 if (isMaster) fRun->EndOfRun(fEdeptrue, fRmstrue, fLimittrue);    

 // show Rndm status
 if (isMaster) G4Random::showEngineStatus();

 // save histos and close analysis
 fAnalysisManager->Write();
 fAnalysisManager->CloseFile(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetEdepAndRMS(G4ThreeVector Value)
{
  fEdeptrue = Value(0);
  fRmstrue  = Value(1);
  fLimittrue= Value(2);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetVerbose(G4int val)  
{
  fVerbose = val;
  if (fRun) fRun->SetVerbose(val);
}
     
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
