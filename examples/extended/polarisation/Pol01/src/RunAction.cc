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
// $Id: RunAction.cc 86443 2014-11-12 09:40:56Z gcosmo $
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
  : detector(det), primary(prim), ProcCounter(0), fAnalysisManager(0)
{
  totalEventCount=0;
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
  
  if (ProcCounter) delete ProcCounter;
  ProcCounter = new ProcessesCount;
  totalEventCount = 0;
  photonStats.Clear();
  electronStats.Clear();
  positronStats.Clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillData(const G4ParticleDefinition* particle,
                         G4double kinEnergy, G4double costheta, 
                         G4double phi,
                         G4double longitudinalPolarization)
{
  G4int id = -1;
  if (particle == fGamma) { 
    photonStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    if(fAnalysisManager) { id = 1; } 
  } else if (particle == fElectron) { 
    electronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
    if(fAnalysisManager) { id = 5; } 
  } else if (particle == fPositron) {
    positronStats.FillData(kinEnergy, costheta, longitudinalPolarization);
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
  size_t nbProc = ProcCounter->size();
  size_t i = 0;
  while ((i<nbProc)&&((*ProcCounter)[i]->GetName()!=procName)) i++;
  if (i == nbProc) ProcCounter->push_back( new OneProcessCount(procName));
  
  (*ProcCounter)[i]->Count();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  G4int  prec = G4cout.precision(5);
    
  G4Material* material = detector->GetMaterial();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = 
                            primary->GetParticleGun()->GetParticleDefinition();
  G4String Particle = particle->GetParticleName();    
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbOfEvents << " "<< Particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(detector->GetBoxSizeZ(),"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //frequency of processes
  G4cout << "\n Process calls frequency --->\n";
  for (size_t i=0; i< ProcCounter->size();i++) {
     G4String procName = (*ProcCounter)[i]->GetName();
     G4int    count    = (*ProcCounter)[i]->GetCounter(); 
     G4cout << "\t" << procName << " = " << count<<"\n";
  }
  
  if (totalEventCount == 0) return;
  
  G4cout<<" Gamma: \n";
  photonStats.PrintResults(totalEventCount);
  G4cout<<" Electron: \n";
  electronStats.PrintResults(totalEventCount);
  G4cout<<" Positron: \n";
  positronStats.PrintResults(totalEventCount);

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
  ++totalEventCount;
  photonStats.EventFinished();
  electronStats.EventFinished();
  positronStats.EventFinished();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::ParticleStatistics()
  : currentNumber(0),
    totalNumber(0), totalNumber2(0),
    sumEnergy(0), sumEnergy2(0),
    sumPolarization(0), sumPolarization2(0),
    sumCosTheta(0), sumCosTheta2(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::ParticleStatistics::~ParticleStatistics()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::EventFinished()
{
  totalNumber+=currentNumber;
  totalNumber2+=currentNumber*currentNumber;
  currentNumber=0;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics:: FillData(G4double kinEnergy, 
                                              G4double costheta,
                                              G4double longitudinalPolarization)
{
  ++currentNumber;
  sumEnergy+=kinEnergy;
  sumEnergy2+=kinEnergy*kinEnergy;
  sumPolarization+=longitudinalPolarization;
  sumPolarization2+=longitudinalPolarization*longitudinalPolarization;
  sumCosTheta+=costheta;
  sumCosTheta2+=costheta*costheta;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::PrintResults(G4int totalNumberOfEvents)
{
  G4cout<<"Mean Number per Event :"
        <<G4double(totalNumber)/G4double(totalNumberOfEvents)<<"\n";
  if (totalNumber==0) totalNumber=1;
  G4double energyMean=sumEnergy/totalNumber;
  G4double energyRms=std::sqrt(sumEnergy2/totalNumber-energyMean*energyMean);
  G4cout<<"Mean Energy :"<< G4BestUnit(energyMean,"Energy")
        <<" +- "<<G4BestUnit(energyRms,"Energy")<<"\n";
  G4double polarizationMean=sumPolarization/totalNumber;
  G4double polarizationRms=
    std::sqrt(sumPolarization2/totalNumber-polarizationMean*polarizationMean);
  G4cout<<"Mean Polarization :"<< polarizationMean
        <<" +- "<<polarizationRms<<"\n";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::ParticleStatistics::Clear()
{
  currentNumber=0;
  totalNumber=totalNumber2=0;
  sumEnergy=sumEnergy2=0;
  sumPolarization=sumPolarization2=0;
  sumCosTheta=sumCosTheta2=0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
