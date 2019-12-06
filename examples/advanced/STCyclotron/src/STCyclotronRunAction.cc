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
//  Author: F. Poignant, floriane.poignant@gmail.com
//
// file STCyclotronRunAction.cc
//
#include "STCyclotronRunAction.hh"
#include "STCyclotronRunActionMessenger.hh"
#include "STCyclotronPrimaryGeneratorAction.hh"
#include "STCyclotronSensitiveTarget.hh"
#include "STCyclotronRun.hh"
#include "STCyclotronDetectorConstruction.hh"
#include "STCyclotronAnalysis.hh"
#include "G4GeneralParticleSource.hh"
#include "G4UserRunAction.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4AccumulableManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <iomanip>
#include <cmath>


STCyclotronRunAction::STCyclotronRunAction(STCyclotronDetectorConstruction*)
{ 

  //----------------------------------------------
  //Set printing event number per each 100 events
  //G4RunManager::GetRunManager()->SetPrintProgress(100);
  //----------------------------------------------
  fMessenger = new STCyclotronRunActionMessenger(this);
  fIrradiationTime = 3.; //in hour
  fIn = 185; //in mm
  fOut = 188; //in mm
  
  //----------------------------------------------
  // Analysis Manager for storage of data in histograms
  //----------------------------------------------

  auto analysisManager = G4AnalysisManager::Instance();
  G4cout << "Using " << analysisManager->GetType() << G4endl;
  analysisManager->SetVerboseLevel(1);
  
  //The edges/bin are default value that may be modified in the file 'Macro/init_parameters.mac'
  //Create TH1D histograms
  analysisManager->CreateH1("H10","Energy of the primary particles when reaching the target (MeV)", 500,12.,19.); //in MeV
  analysisManager->CreateH1("H11","Energy of the primary particle when reaching the foil (MeV)",100,12.,19.); //in MeV
  analysisManager->CreateH1("H12","Energy spectrum of primaries going out from the target (MeV)", 100, 0., 19.); //in MeV
  analysisManager->CreateH1("H13","Energy of the primary particle when going out from the foil (MeV)",100,0,19);
  analysisManager->CreateH1("H14","Depth of isotope creation in the target (mm)", 300, fIn, fOut);
  analysisManager->CreateH1("H15","Energy spectrum of the positrons created in the target by the beam and secondaries (MeV)", 100, 0., 17.); 
  analysisManager->CreateH1("H16","Energy spectrum of the electrons created in the target by the beam and secondaries (MeV)", 100, 0., 17.); //in MeV
  analysisManager->CreateH1("H17","Energy spectrum of the gammas created in the target by the beam and secondaries (MeV)", 100, 0., 17.); //in MeV
  analysisManager->CreateH1("H18","Energy spectrum of the neutrons created in the target by the beam and secondaries (MeV)", 100, 0., 17.); //in MeV
  analysisManager->CreateH1("H19","Energy spectrum of the positrons created in the target by the decay (MeV)", 100, 0., 17.);//, MeV);
  analysisManager->CreateH1("H110","Energy spectrum of the electrons created in the target by the decay (MeV)", 100, 0., 17.);
  analysisManager->CreateH1("H111","Energy spectrum of the gammas created in the target (MeV) by the decay", 100, 0., 17.);
  analysisManager->CreateH1("H112","Energy spectrum of the neutrons created in the target (MeV) by the decay", 100, 0., 17.);
  analysisManager->CreateH1("H113","Energy spectrum of the nu_e created in the target (MeV) by the decay", 100, 0., 17.);
  analysisManager->CreateH1("H114","Energy spectrum of the anti_nu_e created in the target (MeV) by the decay", 100, 0., 17.);

  //Create TH2D histograms
  
  analysisManager->CreateH2("H20", "Beam intensity before hiting the target (mm)",100, -7.5 , 7.5, 100, -7.5, 7.5);
  analysisManager->CreateH2("H21", "Beam intensity before hiting the foil (mm)",100, -7.5 , 7.5, 100, -7.5, 7.5);
  analysisManager->CreateH2("H22", "Radioisotopes produced", 11, 24.5, 35.5, 20, 54.5, 74.5);
  analysisManager->CreateH2("H23", "Energy (MeV) = f(depth (mm))", 100, fIn, fOut,100,0.,16.);
  analysisManager->CreateH2("H24", "Beam intensity going out from the target (mm)",100, -7.5 , 7.5, 100, -7.5, 7.5);
  analysisManager->CreateH2("H25", "Beam intensity going out from the foil (mm)", 100, -7.5 , 7.5, 100, -7.5, 7.5);
}

STCyclotronRunAction::~STCyclotronRunAction()
{
  delete fMessenger;
  delete G4AnalysisManager::Instance();
}

G4Run* STCyclotronRunAction::GenerateRun()
{ 
  fRun = new STCyclotronRun();
  return fRun;
}

void STCyclotronRunAction::BeginOfRunAction(const G4Run*)
{ 
  //----------------------------------------------
  //     Inform the runManager to save random number seed
  //----------------------------------------------
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->OpenFile("SolidTargetCyclotron");
}

void STCyclotronRunAction::EndOfRunAction(const G4Run*)
{
  if(isMaster)fRun->EndOfRun(fIrradiationTime);
  auto analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
 
}

void STCyclotronRunAction::SetIrradiationTime(G4double time)
{
  if(fIrradiationTime != time){
    fIrradiationTime = time;
    G4cout << "The time of irradiation is now the following : " << fIrradiationTime << " hour(s)." << G4endl;
  }
}
