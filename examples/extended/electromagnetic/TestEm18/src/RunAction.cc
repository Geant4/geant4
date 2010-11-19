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
// $Id: RunAction.cc,v 1.3 2010-11-19 12:17:50 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "HistoManager.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4EmCalculator.hh"

#include "Randomize.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin,
                     HistoManager* histo)
:detector(det), primary(kin), histoManager(histo)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* run)
{
  G4cout << "### Run " << run->GetRunID() << " start." << G4endl;

  //initialisation
  //
  energyDeposit = 0.;
  
  nbCharged = nbNeutral = 0;
  energyCharged = energyNeutral = 0.;  
  emin[0] = emin[1] = DBL_MAX;
  emax[0] = emax[1] = 0.;    
    
  nbSteps = 0;
  trackLength = 0.; 

  histoManager->book();

  // do not save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nbEvents = aRun->GetNumberOfEvent();
  if (nbEvents == 0) return;
  
  G4Material* material = detector->GetMaterial();
  G4double length  = detector->GetSize();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double eprimary = primary->GetParticleGun()->GetParticleEnergy();
  
  G4int prec = G4cout.precision(3);
  G4cout << "\n ======================== run summary ======================\n";  
  G4cout << "\n The run was " << nbEvents << " " << partName << " of "
         << G4BestUnit(eprimary,"Energy") << " through " 
	 << G4BestUnit(length,"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")";
  G4cout << "\n ===========================================================\n";
  G4cout << G4endl;
  
  histoManager->save();
   
  G4cout.precision(5);
  
  //track length
  //
  G4double trackLPerEvent = trackLength/nbEvents;
  G4double nbStepPerEvent = G4double(nbSteps)/nbEvents;
  G4double stepSize = trackLength/nbSteps;
  
  G4cout 
    << "\n trackLength= " 
    << G4BestUnit(trackLPerEvent, "Length")
    << "\t nb of steps= " << nbStepPerEvent
    << "  stepSize= " << G4BestUnit(stepSize, "Length")
    << G4endl;
      
  //charged secondaries (ionization, direct pair production)
  //
  G4double energyPerEvent = energyCharged/nbEvents;
  G4double nbPerEvent = G4double(nbCharged)/nbEvents;
  G4double meanEkin = 0.;
  if (nbCharged) meanEkin = energyCharged/nbCharged;
  
  G4cout 
    << "\n d-rays  : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t  nb of d-rays= " << nbPerEvent
    << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
    << "  Tmin= "   << G4BestUnit(emin[0],  "Energy")
    << "  Tmax= "   << G4BestUnit(emax[0],  "Energy")
    << G4endl;
         
  //neutral secondaries (bremsstrahlung)
  //
  energyPerEvent = energyNeutral/nbEvents;
  nbPerEvent = G4double(nbNeutral)/nbEvents;
  meanEkin = 0.;
  if (nbNeutral) meanEkin = energyNeutral/nbNeutral;
  
  G4cout 
    << "\n brems   : eLoss/primary= " 
    << G4BestUnit(energyPerEvent, "Energy")
    << "\t  nb of gammas= " << nbPerEvent
    << "  <Tkin>= " << G4BestUnit(meanEkin, "Energy")
    << "  Tmin= "   << G4BestUnit(emin[1],  "Energy")
    << "  Tmax= "   << G4BestUnit(emax[1],  "Energy")
    << G4endl;
    
  // Computations below only for charged particles
  if (particle->GetPDGCharge() == 0.) return;

  G4EmCalculator emCal;
  
  //local energy deposit
  //
  energyPerEvent = energyDeposit/nbEvents;
  //
  G4double r0  = emCal.GetRangeFromRestricteDEDX(eprimary,particle,material);  
  G4double r1 = r0 - trackLPerEvent;
  G4double etry = eprimary - energyPerEvent;  
  G4double efinal = 0.;
  if (r1 > 0. && etry > 0.0) efinal = GetEnergyFromRestrictedRange(r1,particle,material,etry);
  G4double dEtable = eprimary - efinal;
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
  G4double energyTotal = energyDeposit + energyCharged + energyNeutral;
  energyPerEvent = energyTotal/nbEvents;
  //
  r0  = emCal.GetCSDARange(eprimary,particle,material);  
  r1 = r0 - trackLPerEvent;
  etry = eprimary - energyPerEvent;
  efinal = 0.;
  //G4cout << "r0= " << r0 << "  r1= " << r1 << "  " << particle->GetParticleName()
  //	 << " etry= " << etry << "  " << material->GetName() << " e0= " << eprimary << G4endl;
  if (r1 > 0.0 && etry > 0.0) efinal = GetEnergyFromCSDARange(r1,particle,material,etry);
  dEtable = eprimary - efinal;
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
