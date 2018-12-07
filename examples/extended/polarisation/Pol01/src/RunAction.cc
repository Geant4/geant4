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
/// \file polarisation/Pol01/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4ParticleDefinition.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
: G4UserRunAction(),
  fGamma(0), fElectron(0), fPositron(0),
  fDetector(det), fPrimary(prim), fProcCounter(0), fAnalysisManager(0), 
  fTotalEventCount(0),
  fPhotonStats(), fElectronStats(), fPositronStats()
{
  fGamma = G4Gamma::Gamma();
  fElectron = G4Electron::Electron();
  fPositron = G4Positron::Positron();

  BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  //  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  //  CLHEP::HepRandom::showEngineStatus();
  
  if (fProcCounter) delete fProcCounter;
  fProcCounter = new ProcessesCount;
  fTotalEventCount = 0;
  fPhotonStats.Clear();
  fElectronStats.Clear();
  fPositronStats.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillData(const G4ParticleDefinition* particle,
                         G4double kinEnergy, G4double costheta, 
                         G4double phi,
                         G4double longitudinalPolarization)
{
  G4int id = -1;
  if (particle == fGamma) { 
    fPhotonStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    if(fAnalysisManager) { id = 1; } 
  } else if (particle == fElectron) { 
    fElectronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    if(fAnalysisManager) { id = 5; } 
  } else if (particle == fPositron) {
    fPositronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    if(fAnalysisManager) { id = 9; } 
  }
  if(id > 0) {
    fAnalysisManager->FillH1(id,kinEnergy,1.0);
    fAnalysisManager->FillH1(id+1,costheta,1.0);
    fAnalysisManager->FillH1(id+2,phi,1.0);
    fAnalysisManager->FillH1(id+3,longitudinalPolarization,1.0);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BookHisto()
{
  // Always creating analysis manager
  fAnalysisManager = G4AnalysisManager::Instance();
  fAnalysisManager->SetActivation(true);
  fAnalysisManager->SetVerboseLevel(1);

  // Open file histogram file
  fAnalysisManager->OpenFile("pol01");
  fAnalysisManager->SetFirstHistoId(1);

  // Creating an 1-dimensional histograms in the root directory of the tree
  const G4String id[] = { "h1", "h2", "h3", "h4", "h5", 
                          "h6", "h7", "h8", "h9", "h10", "h11", "h12"};
  const G4String title[] = 
                { "Gamma Energy distribution",                      //1
                  "Gamma Cos(Theta) distribution",                  //2
                  "Gamma Phi angular distribution",                 //3
                  "Gamma longitudinal Polarization",                //4
                  "Electron Energy distribution",                   //5
                  "Electron Cos(Theta) distribution",               //6
                  "Electron Phi angular distribution",              //7
                  "Electron longitudinal Polarization",             //8
                  "Positron Energy distribution",                   //9
                  "Positron Cos(Theta) distribution",               //10
                  "Positron Phi angular distribution",              //11
                  "Positron longitudinal Polarization"              //12
                 };
  G4double vmin, vmax;
  G4int nbins = 120;
  for(int i=0; i<12; ++i) {
    G4int j = i - i/4*4;
    if(0==j)      { vmin = 0.; vmax = 12.*MeV; }
    else if(1==j) { vmin = -1.; vmax = 1.; }
    else if(2==j) { vmin = 0.; vmax = pi; }
    else          { vmin = -1.5; vmax = 1.5; }
    G4int ih = fAnalysisManager->CreateH1(id[i],title[i],nbins,vmin,vmax);
    fAnalysisManager->SetH1Activation(ih, false);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SaveHisto(G4int nevents)
{
  if(fAnalysisManager) {
    G4double norm = 1.0/G4double(nevents);
    for(int i=0; i<12; ++i) {
      fAnalysisManager->ScaleH1(i, norm);
    }
    fAnalysisManager->Write();
    fAnalysisManager->CloseFile();
    delete fAnalysisManager;
    fAnalysisManager = 0;
  } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName)
{
  // is the process already counted ?
  // *AS* change to std::map?!
  size_t nbProc = fProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*fProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) fProcCounter->push_back( new ProcessCount(procName));
  
  (*fProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4int  prec = G4cout.precision(5);
    
  G4Material* material = fDetector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = 
                            fPrimary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(fDetector->GetBoxSizeZ(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  G4cout << "\n Process calls frequency --->\n";
  for (size_t i=0; i< fProcCounter->size();i++) {
     G4String procName = (*fProcCounter)[i]->GetName();
     G4int    count    = (*fProcCounter)[i]->GetCounter(); 
     G4cout << "\t" << procName << " = " << count<<"\n";
  }
  
  if (fTotalEventCount == 0) return;
  
  G4cout<<" Gamma: \n";
  fPhotonStats.PrintResults(fTotalEventCount);
  G4cout<<" Electron: \n";
  fElectronStats.PrintResults(fTotalEventCount);
  G4cout<<" Positron: \n";
  fPositronStats.PrintResults(fTotalEventCount);

  //cross check from G4EmCalculator
  //  G4cout << "\n Verification from G4EmCalculator. \n";  
  //  G4EmCalculator emCal;

  //restore default format         
  G4cout.precision(prec);         

  // write out histograms  
  SaveHisto(NbOfEvents);

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EventFinished()
{
  ++fTotalEventCount;
  fPhotonStats.EventFinished();
  fElectronStats.EventFinished();
  fPositronStats.EventFinished();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::ParticleStatistics()
  : fCurrentNumber(0),
    fTotalNumber(0), fTotalNumber2(0),
    fSumEnergy(0), fSumEnergy2(0),
    fSumPolarization(0), fSumPolarization2(0),
    fSumCosTheta(0), fSumCosTheta2(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::~ParticleStatistics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::EventFinished()
{
  fTotalNumber+=fCurrentNumber;
  fTotalNumber2+=fCurrentNumber*fCurrentNumber;
  fCurrentNumber=0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics:: FillData(G4double kinEnergy, 
                                              G4double costheta,
                                              G4double longitudinalPolarization)
{
  ++fCurrentNumber;
  fSumEnergy+=kinEnergy;
  fSumEnergy2+=kinEnergy*kinEnergy;
  fSumPolarization+=longitudinalPolarization;
  fSumPolarization2+=longitudinalPolarization*longitudinalPolarization;
  fSumCosTheta+=costheta;
  fSumCosTheta2+=costheta*costheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::PrintResults(G4int totalNumberOfEvents)
{
  G4cout<<"Mean Number per Event :"
        <<G4double(fTotalNumber)/G4double(totalNumberOfEvents)<<"\n";
  if (fTotalNumber==0) fTotalNumber=1;
  G4double energyMean=fSumEnergy/fTotalNumber;
  G4double energyRms=std::sqrt(fSumEnergy2/fTotalNumber-energyMean*energyMean);
  G4cout<<"Mean Energy :"<< G4BestUnit(energyMean,"Energy")
        <<" +- "<<G4BestUnit(energyRms,"Energy")<<"\n";
  G4double polarizationMean=fSumPolarization/fTotalNumber;
  G4double polarizationRms=
    std::sqrt(fSumPolarization2/fTotalNumber-polarizationMean*polarizationMean);
  G4cout<<"Mean Polarization :"<< polarizationMean
        <<" +- "<<polarizationRms<<"\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::Clear()
{
  fCurrentNumber=0;
  fTotalNumber=fTotalNumber2=0;
  fSumEnergy=fSumEnergy2=0;
  fSumPolarization=fSumPolarization2=0;
  fSumCosTheta=fSumCosTheta2=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
