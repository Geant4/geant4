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
/// \file electromagnetic/TestEm18/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:G4UserRunAction(),fDetector(det), fPrimary(kin), fHistoManager(0)
{ 
  fHistoManager = new HistoManager(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ 
  delete fHistoManager; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run*)
{
  //initialisation
  //
  fNbSteps = 0;
  fTrackLength = 0.;
  fStepMin = DBL_MAX;
  fStepMax = 0.;

  fEdepPrimary = fEdepSecondary = fEdepTotal = 0.;
  fEdepPrimMin = fEdepSecMin = fEdepTotMin = DBL_MAX;
  fEdepPrimMax = fEdepSecMax = fEdepTotMax = 0.;

  fEnergyTransfered = 0.;
  fEtransfMin = DBL_MAX;
  fEtransfMax = 0.;

  fEnergyLost = 0.;
  fElostMin = DBL_MAX;
  fElostMax = 0.;

  fEnergyBalance = 0.;
  fEbalMin = DBL_MAX;
  fEbalMax = 0.;

  //histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if ( analysisManager->IsActive() ) {
    analysisManager->OpenFile();
  }

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CountProcesses(G4String procName) 
{
  std::map<G4String,G4int>::iterator it = fProcCounter.find(procName);
  if ( it == fProcCounter.end()) {
    fProcCounter[procName] = 1;
  }
  else {
    fProcCounter[procName]++; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TrackLength (G4double step)
{
  fTrackLength += step; fNbSteps++;
  if (step<fStepMin) fStepMin = step;
  if (step>fStepMax) fStepMax = step;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EnergyDeposited (G4double edepPrim, G4double edepSecond)
{
  fEdepPrimary += edepPrim;
  if (edepPrim<fEdepPrimMin) fEdepPrimMin = edepPrim;
  if (edepPrim>fEdepPrimMax) fEdepPrimMax = edepPrim;
  
  fEdepSecondary += edepSecond;
  if (edepSecond<fEdepSecMin) fEdepSecMin = edepSecond;
  if (edepSecond>fEdepSecMax) fEdepSecMax = edepSecond;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EnergyTransferedByProcess(G4String process, G4double energy)
{
  std::map<G4String, MinMaxData>::iterator it = fEtransfByProcess.find(process);
  if ( it == fEtransfByProcess.end()) {
    fEtransfByProcess[process] = MinMaxData(1, energy, energy, energy);
  }
  else {
    MinMaxData& data = it->second;
    data.fCount++;
    data.fVsum += energy;
    //update min max
    G4double emin = data.fVmin;
    if (energy < emin) data.fVmin = energy;
    G4double emax = data.fVmax;
    if (energy > emax) data.fVmax = energy; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EnergyTransfered (G4double energy)
{
  fEnergyTransfered += energy;
  if (energy<fEtransfMin) fEtransfMin = energy;
  if (energy>fEtransfMax) fEtransfMax = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TotalEnergyLost (G4double energy)
{
  fEnergyLost += energy;
  if (energy<fElostMin) fElostMin = energy;
  if (energy>fElostMax) fElostMax = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EnergyBalance (G4double energy)
{
  fEnergyBalance += energy;
  if (energy<fEbalMin) fEbalMin = energy;
  if (energy>fEbalMax) fEbalMax = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::TotalEnergyDeposit (G4double energy)
{
  fEdepTotal += energy;
  if (energy<fEdepTotMin) fEdepTotMin = energy;
  if (energy>fEdepTotMax) fEdepTotMax = energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EnergySpectrumOfSecondaries(G4String particle, G4double energy)
{
 std::map<G4String,MinMaxData>::iterator it = fEkinOfSecondaries.find(particle);
  if ( it == fEkinOfSecondaries.end()) {
    fEkinOfSecondaries[particle] = MinMaxData(1, energy, energy, energy);
  }
  else {
    MinMaxData& data = it->second;
    data.fCount++;
    data.fVsum += energy;
    //update min max
    G4double emin = data.fVmin;
    if (energy < emin) data.fVmin = energy;
    G4double emax = data.fVmax;
    if (energy > emax) data.fVmax = energy; 
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nbEvents = aRun->GetNumberOfEvent();
  if (nbEvents == 0) return;
  
  G4Material* material = fDetector->GetMaterial();
  G4double length  = fDetector->GetSize();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = fPrimary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double ePrimary = fPrimary->GetParticleGun()->GetParticleEnergy();
  
  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";
  G4cout << "\n The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(ePrimary,"Energy") << " through " 
         << G4BestUnit(length,"Length") << " of "
         << material->GetName() << " (density: " 
         << G4BestUnit(density,"Volumic Mass") << ")";
  G4cout << G4endl;

  if (particle->GetPDGCharge() == 0.) return;

  G4cout.precision(4);

  //frequency of processes
  //
  G4cout << "\n Process defining step :" << G4endl;
  G4int index = 0;
  for ( const auto& procCounter : fProcCounter ) {
     G4String procName = procCounter.first;
     G4int    count    = procCounter.second;
     G4String space = " "; if (++index%4 == 0) space = "\n";
     G4cout << " " << std::setw(15) << procName << "="<< std::setw(7) << count
            << space;
  }
  G4cout << G4endl;

  //track length
  //
  G4double trackLPerEvent = fTrackLength/nbEvents;
  G4double nbStepPerEvent = double(fNbSteps)/nbEvents;
  G4double stepSize = fTrackLength/fNbSteps;
  
  G4cout 
    << "\n TrackLength = " 
    << G4BestUnit(trackLPerEvent, "Length")
    << "  nb of steps = " << nbStepPerEvent
    << "  stepSize = " << G4BestUnit(stepSize, "Length")
    << "  ("  << G4BestUnit(fStepMin, "Length")
    << "--> " << G4BestUnit(fStepMax, "Length") << ")"
    << G4endl;

  //continuous energy deposited by primary track dE1
  //
  G4double energyPerEvent = fEdepPrimary/nbEvents;
  
  G4cout 
    << "\n Energy continuously deposited along primary track"
    << " (restricted dE/dx)  dE1 = "
    << G4BestUnit(energyPerEvent, "Energy")
    << "  ("   << G4BestUnit(fEdepPrimMin, "Energy")
    << " --> " << G4BestUnit(fEdepPrimMax, "Energy") << ")"
    << G4endl;
  
  //eveluation of dE1 from reading restricted Range table
  //
  G4EmCalculator emCal;
  
  G4double r0  = emCal.GetRangeFromRestricteDEDX(ePrimary,particle,material);
  G4double r1 = r0 - trackLPerEvent;
  G4double etry = ePrimary - energyPerEvent;  
  G4double efinal = 0.;
  if (r1 > 0.) efinal = GetEnergyFromRestrictedRange(r1,particle,material,etry);
  G4double dEtable = ePrimary - efinal;
  G4double ratio = 0.;
  if (dEtable > 0.) ratio = energyPerEvent/dEtable;
    
  G4cout 
    << "\n Evaluation of dE1 from reading restricted Range table : dE1_table = "
    << G4BestUnit(dEtable, "Energy")
    << "   ---> dE1/dE1_table = " << ratio
    << G4endl;
  
  // energy transfered to secondary particles by process : dE2
  //
  G4cout << "\n Energy transfered to secondary particles :" << G4endl;
  std::map<G4String,MinMaxData>::iterator it1;               
  for (it1 = fEtransfByProcess.begin(); it1 != fEtransfByProcess.end(); it1++) {
     G4String name = it1->first;
     MinMaxData data = it1->second;
     energyPerEvent = data.fVsum/nbEvents;
     G4double eMin = data.fVmin;
     G4double eMax = data.fVmax;

     G4cout << "  " << std::setw(17) << "due to " + name << ":  dE2 = "
            << std::setw(6) << G4BestUnit(energyPerEvent, "Energy")
            << "  ("   << G4BestUnit(eMin, "Energy")
            << " --> " << G4BestUnit(eMax, "Energy")
            << ")" << G4endl;           
  }

  // total energy tranfered : dE3 = sum of dE2
  //
  energyPerEvent = fEnergyTransfered/nbEvents;
  
  G4cout 
    << "\n Total energy transfered to secondaries : dE3 = sum of dE2 = "
    << G4BestUnit(energyPerEvent, "Energy")
    << "  ("   << G4BestUnit(fEtransfMin, "Energy")
    << " --> " << G4BestUnit(fEtransfMax, "Energy") << ")"
    << G4endl;

  // total energy lost by incident particle : dE4 = dE1 + dE3
  //
  energyPerEvent = fEnergyLost/nbEvents;
  
  G4cout 
    << "\n Total energy lost by incident particle : dE4 = dE1 + dE3 = "
    << G4BestUnit(energyPerEvent, "Energy")
    << "  ("   << G4BestUnit(fElostMin, "Energy")
    << " --> " << G4BestUnit(fElostMax, "Energy") << ")"
    << G4endl;
  
  // calcul of energy lost from energy balance : dE4_bal = E_in - E_out
  //
  energyPerEvent = fEnergyBalance/nbEvents;
  
  G4cout 
    << "\n calcul of dE4 from energy balance : dE4_bal = E_in - E_out = "
    << G4BestUnit(energyPerEvent, "Energy")
    << "  ("   << G4BestUnit(fEbalMin, "Energy")
    << " --> " << G4BestUnit(fEbalMax, "Energy") << ")"
    << G4endl;
  
  //eveluation of dE4 from reading full Range table
  //
  r0  = emCal.GetCSDARange(ePrimary,particle,material);
  r1 = r0 - trackLPerEvent;
  etry = ePrimary - energyPerEvent;  
  efinal = 0.;
  if (r1 > 0.) efinal = GetEnergyFromCSDARange(r1,particle,material,etry);
  dEtable = ePrimary - efinal;
  ratio = 0.;
  if (dEtable > 0.) ratio = energyPerEvent/dEtable;
    
  G4cout 
    << "\n Evaluation of dE4 from reading full Range table : dE4_table = "
    << G4BestUnit(dEtable, "Energy")
    << "   ---> dE4/dE4_table = " << ratio
    << G4endl;
  
  //energy spectrum of secondary particles
  //
  G4cout << "\n Energy spectrum of secondary particles :" << G4endl;
  std::map<G4String,MinMaxData>::iterator it2;               
  for (it2 = fEkinOfSecondaries.begin();it2 != fEkinOfSecondaries.end(); it2++){
     G4String name = it2->first;
     MinMaxData data = it2->second;
     G4int count = data.fCount;
     G4double eMean = data.fVsum/count;
     G4double eMin = data.fVmin;
     G4double eMax = data.fVmax;

     G4cout << "  " << std::setw(13) << name << ": " << std::setw(7) << count
            << "  Emean = " << std::setw(6) << G4BestUnit(eMean, "Energy")
            << "  ("   << G4BestUnit(eMin, "Energy")
            << " --> " << G4BestUnit(eMax, "Energy")
            << ")" << G4endl;           
  }
  G4cout << G4endl;
  
  //continuous energy deposited by secondary tracks dE5
  // (only if secondary particles are tracked)
  //
  if (fEdepSecondary > 0.) {
    energyPerEvent = fEdepSecondary/nbEvents;

    G4cout 
      << "\n Energy continuously deposited along secondary tracks"
      << " (restricted dE/dx)  dE5 = "
      << G4BestUnit(energyPerEvent, "Energy")
      << "  ("   << G4BestUnit(fEdepSecMin, "Energy")
      << " --> " << G4BestUnit(fEdepSecMax, "Energy") << ")"
      << G4endl;

    // total energy deposited : dE6 = dE1 + dE5
    //
    energyPerEvent = fEdepTotal/nbEvents;

    G4cout 
      << "\n Total energy deposited : dE6 = dE1 + dE5 = "
      << G4BestUnit(energyPerEvent, "Energy")
      << "  ("   << G4BestUnit(fEdepTotMin, "Energy")
      << " --> " << G4BestUnit(fEdepTotMax, "Energy") << ") \n"
      << G4endl;
}
  
  G4cout.precision(prec);
  
  //clear maps
  //
  fProcCounter.clear();
  fEtransfByProcess.clear();
  fEkinOfSecondaries.clear();

  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromRestrictedRange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetRangeFromRestricteDEDX(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.GetDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;    
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromRestRange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::GetEnergyFromCSDARange(G4double range,
            G4ParticleDefinition* particle, G4Material* material, G4double Etry)
{
  G4EmCalculator emCal;
    
  G4double Energy = Etry, dE = 0., dEdx;
  G4double r, dr;
  G4double err  = 1., errmax = 0.00001;
  G4int    iter = 0 , itermax = 10;  
  while (err > errmax && iter < itermax) {
    iter++;
    Energy -= dE;
    r = emCal.GetCSDARange(Energy,particle,material);
    dr = r - range;          
    dEdx = emCal.ComputeTotalDEDX(Energy,particle,material);
    dE = dEdx*dr;
    err = std::abs(dE)/Energy;
  }
  if (iter == itermax) {
    G4cout 
    << "\n  ---> warning: RunAction::GetEnergyFromCSDARange() did not converge"
    << "   Etry = " << G4BestUnit(Etry,"Energy")
    << "   Energy = " << G4BestUnit(Energy,"Energy")
    << "   err = " << err
    << "   iter = " << iter << G4endl;
  }         
         
  return Energy;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
