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
/// \file electromagnetic/TestEm17/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"
#include "MuCrossSections.hh"
#include "G4ProductionCutsTable.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim,
                     HistoManager* HistM)
  : G4UserRunAction(),
    fDetector(det), fPrimary(prim), fProcCounter(0), fHistoManager(HistM)
{
  fMucs = new MuCrossSections();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fMucs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  CLHEP::HepRandom::showEngineStatus();

  fProcCounter = new ProcessesCount();
  fHistoManager->Book();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(const G4String& procName)
{
   //does the process  already encounted ?
   size_t n = fProcCounter->size();
   for(size_t i = 0; i<n; ++i) {
     if((*fProcCounter)[i]->GetName()==procName) {
       (*fProcCounter)[i]->Count();
       return;
     }
   }
   OneProcessCount* count = new OneProcessCount(procName);
   count->Count();
   fProcCounter->push_back(count);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  if (NbOfEvents == 0) return;
  
  //  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
    
  const G4Material* material = fDetector->GetMaterial();
  G4double length  = fDetector->GetSize();
  G4double density = material->GetDensity();
   
  G4String particle = fPrimary->GetParticleGun()->GetParticleDefinition()
                      ->GetParticleName();    
  G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
  
  G4cout << "\n The run consists of " << NbOfEvents << " "<< particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  //total number of process calls
  G4double countTot = 0.;
  G4cout << "\n Number of process calls --->";
  for (size_t i=0; i< fProcCounter->size();++i) {
     G4String procName = (*fProcCounter)[i]->GetName();
     if (procName != "Transportation") {
       G4int count = (*fProcCounter)[i]->GetCounter(); 
       G4cout << "\t" << procName << " : " << count;
       countTot += count;
     }
  }
  
  //compute totalCrossSection, meanFreePath and massicCrossSection
  //
  G4double totalCrossSection = countTot/(NbOfEvents*length);
  G4double MeanFreePath      = 1./totalCrossSection;        
  G4double massCrossSection  = totalCrossSection/density;     
   
  G4cout.precision(5);
  G4cout << "\n Simulation: "
         <<    "total CrossSection = " << totalCrossSection*cm << " /cm"
         << "\t MeanFreePath = "       << G4BestUnit(MeanFreePath,"Length")
         << "\t massicCrossSection = " << massCrossSection*g/cm2 << " cm2/g"
         << G4endl;
  
  //compute theoretical predictions
  //
  if(particle == "mu+" || particle == "mu-") { 
    totalCrossSection = 0.;
    for (size_t i=0; i< fProcCounter->size();++i) {
      G4String procName = (*fProcCounter)[i]->GetName();
      if (procName != "Transportation") {
        totalCrossSection += ComputeTheory(procName, NbOfEvents);
        FillCrossSectionHisto(procName, NbOfEvents);
      }
    }
  
    MeanFreePath     = 1./totalCrossSection;
    massCrossSection = totalCrossSection/density;
  
    G4cout << " Theory:     "
           <<    "total CrossSection = " << totalCrossSection*cm << " /cm"
           << "\t MeanFreePath = "       << G4BestUnit(MeanFreePath,"Length")
           << "\t massicCrossSection = " << massCrossSection*g/cm2 << " cm2/g"
           << G4endl;
  }
                                                                            
  //  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);         

  // delete and remove all contents in fProcCounter 
  size_t n = fProcCounter->size();
  for(size_t i = 0; i<n; ++i) { delete (*fProcCounter)[i]; }
  delete fProcCounter;
  
  fHistoManager->Save();
  
  // show Rndm status
  //CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::ComputeTheory(const G4String& process, G4int NbOfMu)    
{   
  const G4Material* material = fDetector->GetMaterial();
  G4double ekin = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4double particleMass = fPrimary->GetParticleGun()->GetParticleDefinition()->GetPDGMass();

  G4int id = 0; G4double cut = 1.e-10*ekin;
  if (process == "muIoni")          {id = 11; cut =  GetEnergyCut(material,1);}
  else if (process == "muPairProd") {id = 12; cut = 2*(GetEnergyCut(material,1)
                                                      + electron_mass_c2); }
  else if (process == "muBrems")    {id = 13; cut =  GetEnergyCut(material,0);}
  else if (process == "muonNuclear"){id = 14; cut = 100*MeV;}
  else if (process == "muToMuonPairProd"){id = 18; cut = 2*particleMass;}
  if (id == 0) { return 0.; }
  
  G4int nbOfBins = 100;
  //G4double binMin = -10.;
  G4double binMin = std::log10(cut/ekin);
  G4double binMax = 0.;
  G4double binWidth = (binMax-binMin)/G4double(nbOfBins);

  //create histo for theoretical crossSections, with same bining as simulation
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  G4H1* histoTh = 0;
  if (fHistoManager->HistoExist(id)) {
    histoTh  = analysisManager->GetH1(fHistoManager->GetHistoID(id));  
    nbOfBins = fHistoManager->GetNbins(id);
    binMin   = fHistoManager->GetVmin (id);
    binMax   = fHistoManager->GetVmax (id);
    binWidth = fHistoManager->GetBinWidth(id);    
  }
  
  //compute and plot differential crossSection, as function of energy transfert.
  //compute and return integrated crossSection for a given process.
  //(note: to compare with simulation, the integrated crossSection is function
  //       of the energy cut.) 
  // 
  G4double lgeps, etransf, sigmaE, dsigma;
  G4double sigmaTot = 0.;
  const G4double ln10 = std::log(10.);  
  G4double length = fDetector->GetSize();

  //G4cout << "MU: " << process << " E= " << ekin 
  //       <<"  binMin= " << binMin << " binW= " << binWidth << G4endl;

  for (G4int ibin=0; ibin<nbOfBins; ibin++) {
    lgeps = binMin + (ibin+0.5)*binWidth;
    etransf = ekin*std::pow(10.,lgeps);
    sigmaE = fMucs->CR_Macroscopic(process,material,ekin,etransf);
    dsigma = sigmaE*etransf*binWidth*ln10;
    if (etransf > cut) sigmaTot += dsigma;    
    if (histoTh) {
      G4double NbProcess = NbOfMu*length*dsigma;
      histoTh->fill(lgeps, NbProcess);
    }
  }
     
  //return integrated crossSection
  //
  return sigmaTot;   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillCrossSectionHisto(const G4String& process, G4int)
{
  const G4Material* material = fDetector->GetMaterial();
  G4double ekin = fPrimary->GetParticleGun()->GetParticleEnergy();
  G4ParticleDefinition *particle = fPrimary->GetParticleGun()->GetParticleDefinition();
  G4double particleMass = particle->GetPDGMass();
  
  G4EmCalculator emCal;

  G4int id = 0; G4double cut = 1.e-10*ekin;
  if (process == "muIoni")          {id = 21; cut = GetEnergyCut(material,1);}
  else if (process == "muPairProd") {id = 22; cut = 2*(GetEnergyCut(material,1)
                                                      + electron_mass_c2); }
  else if (process == "muBrems")    {id = 23; cut = GetEnergyCut(material,0);}
  else if (process == "muonNuclear"){id = 24; cut = 100*MeV;}
  else if (process == "muToMuonPairProd"){id = 28; cut = 2*particleMass;}
  if (id == 0) { return; }

  G4int nbOfBins = 100;
  G4double binMin = cut;
  G4double binMax = ekin;
  G4double binWidth = (binMax-binMin)/G4double(nbOfBins);

  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  G4H1* histoTh = 0;
  if (fHistoManager->HistoExist(id)) {
    histoTh  = analysisManager->GetH1(fHistoManager->GetHistoID(id));  
    nbOfBins = fHistoManager->GetNbins(id);
    binMin   = fHistoManager->GetVmin (id);
    binMax   = fHistoManager->GetVmax (id);
    binWidth = fHistoManager->GetBinWidth(id);    
  }

  G4double sigma, primaryEnergy;

  for(G4int ibin=0; ibin<nbOfBins; ibin++){
    primaryEnergy = binMin + (ibin+0.5)*binWidth;
    sigma = emCal.GetCrossSectionPerVolume(primaryEnergy, particle, process, material);
    if (histoTh) {
      histoTh->fill(primaryEnergy, sigma);
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyCut(const G4Material* material, G4int idParticle)
{ 
 G4ProductionCutsTable* table = G4ProductionCutsTable::GetProductionCutsTable();
 
 size_t index = 0;
 while ( (table->GetMaterialCutsCouple(index)->GetMaterial() != material) &&
        (index < table->GetTableSize())) index++;

 return (*(table->GetEnergyCutsVector(idParticle)))[index];
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
                   
