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
<<<<<<< HEAD
// $Id: RunAction.cc 68036 2013-03-13 14:13:45Z gcosmo $
=======
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
<<<<<<< HEAD
=======
#include "Randomize.hh"
#include "G4ProductionCutsTable.hh"
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c

#include "Histo.hh"
#include "G4EmCalculator.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
<<<<<<< HEAD
  :G4UserRunAction(),
   fHisto(0),fDetector(det),fKinematic(kin),fRunActionMessenger(0)
{ 
  fVerboseLevel = 1;
  fProjRange = fProjRange2 = fBinLength = fOffsetX = 0.;
  fHisto = new Histo();
  fHisto->SetFileName("monopole");
  // create commands for interactive definition of the detector  
  fRunActionMessenger = new RunActionMessenger(this);
=======
  :fDetector(det),fKinematic(kin)
{
  fMessenger = new RunActionMessenger(this);
  fBinLength = 5 * CLHEP::mm;
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  analysisManager->SetFileName("monopole");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetActivation(true);
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
<<<<<<< HEAD
  delete fHisto;
  delete fRunActionMessenger;
=======
  if(isMaster && G4Threading::IsMultithreadedApplication()) delete fKinematic;

  delete fMessenger;
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  if(GetVerbose() > 0) {
    G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  }
     
  //initialize projected range, tallies, Ebeam, and book histograms
  fProjRange = fProjRange2 = 0.;

  G4double length  = fDetector->GetAbsorSizeX();
  if(0.0 == fBinLength) { fBinLength = 5 * mm; }
  if(fBinLength > fDetector->GetMaxStepSize()) { 
    fBinLength = fDetector->GetMaxStepSize();
  }
  fOffsetX = -0.5 * length;
 
  G4int nbBins = G4lrint(length / fBinLength);

<<<<<<< HEAD

  // Create histograms
  fHisto->Add1D("1","Edep (MeV/mm) along absorber (mm)", nbBins, 0, length, mm);
  fHisto->Add1D("2","DEDX (MeV/mm) of proton", 100, -3., 7.);
  fHisto->Add1D("3","DEDX (MeV/mm) of monopole", 100, -3., 7.);
  fHisto->Add1D("4","Range(mm) of proton", 100, -3., 7., mm);
  fHisto->Add1D("5","Range(mm) of monopole", 100, -3., 7., mm);

  fHisto->Book();
=======
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // Dump production cuts
  G4ProductionCutsTable::GetProductionCutsTable()->DumpCouples();

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  //histograms
  //        
  Book();
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nEvents = aRun->GetNumberOfEvent();
  if (nEvents == 0) { return; }

  //run conditions
  //  
  const G4Material* material = fDetector->GetAbsorMaterial();
  G4double density = material->GetDensity();
  G4String matName = material->GetName();
  const G4ParticleDefinition* part = 
    fKinematic->GetParticleGun()->GetParticleDefinition();
  G4String particle = part->GetParticleName();    
  G4double energy = fKinematic->GetParticleGun()->GetParticleEnergy();

  if(GetVerbose() > 0){
    G4cout << "\n The run consists of " << nEvents << " "<< particle << " of "
           << G4BestUnit(energy,"Energy") << " through " 
           << G4BestUnit(fDetector->GetAbsorSizeX(),"Length") << " of "
           << matName << " (density: " 
           << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  };
         
  //compute projected range and straggling

  fProjRange /= nEvents; fProjRange2 /= nEvents;
  G4double rms = fProjRange2 - fProjRange*fProjRange;        
  if (rms>0.) { rms = std::sqrt(rms); } 
  else { rms = 0.; }

  if(GetVerbose() > 0){
    G4cout.precision(5);       
    G4cout << "\n projected Range= " << G4BestUnit(fProjRange, "Length")
           << "   rms= "             << G4BestUnit(rms, "Length")
           << G4endl;
  };

  G4double ekin[100], dedxproton[100], dedxmp[100];
  G4EmCalculator calc;
  calc.SetVerbose(0);
  G4int i;
  for(i = 0; i < 100; ++i) {
    ekin[i] = std::pow(10., 0.1*G4double(i)) * keV;
    dedxproton[i] = 
      calc.ComputeElectronicDEDX(ekin[i], "proton", matName);
    dedxmp[i] = 
      calc.ComputeElectronicDEDX(ekin[i], "monopole", matName);
  }

  if(GetVerbose() > 0){
    G4cout << "### Stopping Powers" << G4endl;
    for(i=0; i<100; i++) {
      G4cout << " E(MeV)= " << ekin[i] << "  dedxp= " << dedxproton[i]
             << " dedxmp= " << dedxmp[i]
             << G4endl;
    }
  }
  G4cout << "### End of stopping power table" << G4endl;

  // normalize histogram
  G4double fac = (mm/MeV) / (nEvents * fBinLength);
  fHisto->ScaleH1(0,fac);

  if(GetVerbose() > 0){
    G4cout << "Range table for " << matName << G4endl;
  }

  for(i=0; i<100; ++i) {
    G4double e = std::log10(ekin[i] / MeV) + 0.05;
    fHisto->Fill(1, e, dedxproton[i]);
    fHisto->Fill(2, e, dedxmp[i]);
    fHisto->Fill(3, e, std::log10(calc.GetRange(ekin[i],"proton",matName)/mm));
    fHisto->Fill(4, e, std::log10(calc.GetRange(ekin[i],"monopole",matName)/mm));
  }
  
  // save and clean histo
  fHisto->Save();
 
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillHisto(G4int ih, G4double x, G4double weight)
{
<<<<<<< HEAD
  if(GetVerbose() > 1) {
    G4cout << "FillHisto " << ih << "  x=" << x << " weight= " << weight 
           << G4endl;
  }
  fHisto->Fill(ih, x, weight);
=======
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  analysisManager->SetFirstHistoId(1);   

  G4double length = fDetector->GetAbsorSizeX();
  G4int nbBins = G4lrint(length / fBinLength);

  // Create histograms
  analysisManager->CreateH1("h1","Edep (MeV/mm) along absorber (mm)", 
                            nbBins, 0, length);
  analysisManager->CreateH1("h2","Total DEDX (MeV/mm) of proton",100,-3.,7.);
  analysisManager->CreateH1("h3","Total DEDX (MeV/mm) of monopole",100,-3., 7.);
  analysisManager->CreateH1("h4","Range(mm) of proton", 100, -3., 7., "mm");
  analysisManager->CreateH1("h5","Range(mm) of monopole", 100, -3., 7., "mm");
  analysisManager->CreateH1("h6","Restricted DEDX (MeV/mm) of proton",
                            100,-3.,7.);
  analysisManager->CreateH1("h7","Restricted DEDX (MeV/mm) of monopole",
                            100,-3., 7.);
  analysisManager->CreateH1("h8","Delta-electron x-section (1/mm) of proton", 
                            100, -3., 7., "mm");
  analysisManager->CreateH1("h9","Delta-electron x-section (1/mm) of monopole", 
                            100, -3., 7., "mm");
  analysisManager->OpenFile(); 
>>>>>>> 5baee230e93612916bcea11ebf822756cfa7282c
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
