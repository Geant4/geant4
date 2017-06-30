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
// $Id: RunAction.cc 103621 2017-04-19 13:21:45Z gcosmo $
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
  fEnergyDeposit = 0.;
  
  fNbCharged = fNbNeutral = 0;
  fEnergyCharged = fEnergyNeutral = 0.;  
  fEmin[0] = fEmin[1] = DBL_MAX;
  fEmax[0] = fEmax[1] = 0.;    
    
  fNbSteps = 0;
  fTrackLength = 0.;
   
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
  G4cout << "\n ===========================================================\n";
  G4cout << G4endl;
  
  //save histograms      
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  
  if ( analysisManager->IsActive() ) {
    analysisManager->Write();
    analysisManager->CloseFile();
  }      
    
  if (particle->GetPDGCharge() == 0.) return;
   
  G4cout.precision(5);
  
  //track length
  //
  G4double trackLPerEvent = fTrackLength/nbEvents;
  G4double nbStepPerEvent = double(fNbSteps)/nbEvents;
  G4double stepSize = fTrackLength/fNbSteps;
  
  G4cout 
    << "\n TrackLength= " 
    << G4BestUnit(trackLPerEvent, "Length")
    << "\t nb of steps= " << nbStepPerEvent
    << "  stepSize= " << G4BestUnit(stepSize, "Length")
    << G4endl;
      
  //charged secondaries (ionization, direct pair production)
  //
  G4double energyPerEvent = fEnergyCharged/nbEvents;
  G4double nbPerEvent = double(fNbCharged)/nbEvents;
  G4double meanEkin = 0.;
  if (fNbCharged) meanEkin = fEnergyCharged/fNbCharged;
  
  G4cout 
    << "\n d-rays  : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t  nb of d-rays= " << nbPerEvent
    << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
    << "  Tmin= "   << G4BestUnit(fEmin[0], "Energy")
    << "  Tmax= "   << G4BestUnit(fEmax[0], "Energy")
    << G4endl;
         
  //neutral secondaries (bremsstrahlung, pixe)
  //
  energyPerEvent = fEnergyNeutral/nbEvents;
  nbPerEvent = double(fNbNeutral)/nbEvents;
  meanEkin = 0.;
  if (fNbNeutral) meanEkin = fEnergyNeutral/fNbNeutral;
  
  G4cout 
    << "\n gamma   : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t  nb of gammas= " << nbPerEvent
    << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
    << "  Tmin= "   << G4BestUnit(fEmin[1],  "Energy")
    << "  Tmax= "   << G4BestUnit(fEmax[1],  "Energy")
    << G4endl;
    

  G4EmCalculator emCal;
  
  //local energy deposit
  //
  energyPerEvent = fEnergyDeposit/nbEvents;
  //
  G4double r0  = emCal.GetRangeFromRestricteDEDX(ePrimary,particle,material);  
  G4double r1 = r0 - trackLPerEvent;
  G4double etry = ePrimary - energyPerEvent;  
  G4double efinal = 0.;
  if (r1 > 0.) efinal = GetEnergyFromRestrictedRange(r1,particle,material,etry);
  G4double dEtable = ePrimary - efinal;
  G4double ratio = 0.;
  if (dEtable > 0.) ratio = energyPerEvent/dEtable;
    
  G4cout 
    << "\n deposit : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t <dEcut > table= " 
    << G4BestUnit(dEtable, "Energy")
    << "   ---> simul/reference= " << ratio           
    << G4endl;
    
  //total energy transferred
  //
  G4double energyTotal = fEnergyDeposit + fEnergyCharged + fEnergyNeutral;
  energyPerEvent = energyTotal/nbEvents;
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
    << "\n total   : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t <dEfull> table= " 
    << G4BestUnit(dEtable, "Energy")
    << "   ---> simul/reference= " << ratio           
    << G4endl; 

  G4cout.precision(prec);

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
