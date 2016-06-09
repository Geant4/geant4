//
// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
// $Id: RunAction.cc,v 1.20 2005/12/06 11:46:46 gcosmo Exp $
// GEANT4 tag $Name: geant4-08-00 $
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

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //initialisation
  EnergyDeposit  = EnergyDeposit2  = 0.;
  TrakLenCharged = TrakLenCharged2 = 0.;
  TrakLenNeutral = TrakLenNeutral2 = 0.;
  nbStepsCharged = nbStepsCharged2 = 0.;
  nbStepsNeutral = nbStepsNeutral2 = 0.;
  MscProjecTheta = MscProjecTheta2 = 0.;

  nbGamma = nbElect = nbPosit = 0;

  Transmit[0] = Transmit[1] = Reflect[0] = Reflect[1] = 0;

  histoManager->book();

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  // compute mean and rms
  //
  G4int TotNbofEvents = aRun->GetNumberOfEvent();
  if (TotNbofEvents == 0) return;

  EnergyDeposit /= TotNbofEvents; EnergyDeposit2 /= TotNbofEvents;
  G4double rmsEdep = EnergyDeposit2 - EnergyDeposit*EnergyDeposit;
  if (rmsEdep>0.) rmsEdep = std::sqrt(rmsEdep/TotNbofEvents);
  else            rmsEdep = 0.;

  TrakLenCharged /= TotNbofEvents; TrakLenCharged2 /= TotNbofEvents;
  G4double rmsTLCh = TrakLenCharged2 - TrakLenCharged*TrakLenCharged;
  if (rmsTLCh>0.) rmsTLCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsTLCh = 0.;

  TrakLenNeutral /= TotNbofEvents; TrakLenNeutral2 /= TotNbofEvents;
  G4double rmsTLNe = TrakLenNeutral2 - TrakLenNeutral*TrakLenNeutral;
  if (rmsTLNe>0.) rmsTLNe = std::sqrt(rmsTLNe/TotNbofEvents);
  else            rmsTLNe = 0.;

  nbStepsCharged /= TotNbofEvents; nbStepsCharged2 /= TotNbofEvents;
  G4double rmsStCh = nbStepsCharged2 - nbStepsCharged*nbStepsCharged;
  if (rmsStCh>0.) rmsStCh = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStCh = 0.;

  nbStepsNeutral /= TotNbofEvents; nbStepsNeutral2 /= TotNbofEvents;
  G4double rmsStNe = nbStepsNeutral2 - nbStepsNeutral*nbStepsNeutral;
  if (rmsStNe>0.) rmsStNe = std::sqrt(rmsTLCh/TotNbofEvents);
  else            rmsStNe = 0.;

  G4double Gamma = (double)nbGamma/TotNbofEvents;
  G4double Elect = (double)nbElect/TotNbofEvents;
  G4double Posit = (double)nbPosit/TotNbofEvents;

  G4double transmit[2];
  transmit[0] = 100.*Transmit[0]/TotNbofEvents;
  transmit[1] = 100.*Transmit[1]/TotNbofEvents;

  G4double reflect[2];
  reflect[0] = 100.*Reflect[0]/TotNbofEvents;
  reflect[1] = 100.*Reflect[1]/TotNbofEvents;

  G4double rmsMsc = 0.;
  if (Transmit[1] > 0) {
    MscProjecTheta /= (2*Transmit[1]); MscProjecTheta2 /= (2*Transmit[1]);
    rmsMsc = MscProjecTheta2 - MscProjecTheta*MscProjecTheta;
    if (rmsMsc > 0.) rmsMsc = std::sqrt(rmsMsc);
  }
  
  //Stopping Power from input Table.
  //
  G4Material* material = detector->GetAbsorberMaterial();
  G4double length  = detector->GetAbsorberThickness();
  G4double density = material->GetDensity();
   
  G4ParticleDefinition* particle = primary->GetParticleGun()
                                          ->GetParticleDefinition();
  G4String partName = particle->GetParticleName();
  G4double energy = primary->GetParticleGun()->GetParticleEnergy();

  G4EmCalculator emCalculator;
  G4double dEdxTable = 0.;
  if (particle->GetPDGCharge()!= 0.) { 
    dEdxTable = emCalculator.GetDEDX(energy,particle,material);
  }
  G4double stopTable = dEdxTable/density;
  
  //Stopping Power from simulation.
  //    
  G4double meandEdx  = EnergyDeposit/length;
  G4double stopPower = meandEdx/density;  

  G4cout << "\n ======================== run summary ======================\n";

  G4int prec = G4cout.precision(2);
  
  G4cout << "\n The run was " << TotNbofEvents << " " << partName << " of "
         << G4BestUnit(energy,"Energy") << " through " 
	 << G4BestUnit(length,"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  
  G4cout.precision(4);
  G4cout << "\n Total energy deposit in absorber per event = "
         << G4BestUnit(EnergyDeposit,"Energy") << " +- "
         << G4BestUnit(rmsEdep,      "Energy") 
         << G4endl;
	 
  G4cout << "\n Mean dE/dx  = " << meandEdx/(MeV/cm) << " MeV/cm"
         << "\t stopping Power = " << stopPower/(MeV*cm2/g) << " MeV*cm2/g"
	 << G4endl;
	 	 
  G4cout << " (from Table = " << dEdxTable/(MeV/cm) << " MeV/cm)"
         << "\t (from Table    = " << stopTable/(MeV*cm2/g) << " MeV*cm2/g)"
	 << G4endl;
	 
  G4cout << "\n Total track length (charged) in absorber per event = "
         << G4BestUnit(TrakLenCharged,"Length") << " +- "
         << G4BestUnit(rmsTLCh,       "Length") << G4endl;

  G4cout << " Total track length (neutral) in absorber per event = "
         << G4BestUnit(TrakLenNeutral,"Length") << " +- "
         << G4BestUnit(rmsTLNe,       "Length") << G4endl;

  G4cout << "\n Number of steps (charged) in absorber per event = "
         << nbStepsCharged << " +- " << rmsStCh << G4endl;

  G4cout << " Number of steps (neutral) in absorber per event = "
         << nbStepsNeutral << " +- " << rmsStNe << G4endl;

  G4cout << "\n Number of secondaries per event : Gammas = " << Gamma
         << ";   electrons = " << Elect
  	 << ";   positrons = " << Posit << G4endl;

  G4cout << "\n Number of events with the primary particle transmitted = "
         << transmit[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle transmitted "
         << "(same charge as primary) = " << transmit[0] << " %" << G4endl;

  G4cout << "\n Number of events with the primary particle reflected = "
         << reflect[1] << " %" << G4endl;

  G4cout << " Number of events with at least  1 particle reflected "
         << "(same charge as primary) = " << reflect[0] << " %" << G4endl;

  // compute width of the Gaussian central part of the MultipleScattering
  //
  if (histoManager->HistoExist(6)) {
    G4cout << "\n MultipleScattering: rms proj angle of transmit primary particle = "
         << rmsMsc/mrad << " mrad " << G4endl;

    G4cout << " MultipleScattering: computed theta0 (Highland formula)          = "
 	 << ComputeMscHighland()/mrad << " mrad" << G4endl;
  }	 

  G4cout.precision(prec);

  histoManager->save();

  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double RunAction::ComputeMscHighland()
{
  //compute the width of the Gaussian central part of the MultipleScattering
  //projected angular distribution.
  //Eur. Phys. Jour. C15 (2000) page 166, formule 23.9

  G4double t = (detector->GetAbsorberThickness())
              /(detector->GetAbsorberMaterial()->GetRadlen());
  if (t < DBL_MIN) return 0.;

  G4ParticleGun* particle = primary->GetParticleGun();
  G4double T = particle->GetParticleEnergy();
  G4double M = particle->GetParticleDefinition()->GetPDGMass();
  G4double z = std::abs(particle->GetParticleDefinition()->GetPDGCharge()/eplus);

  G4double bpc = T*(T+2*M)/(T+M);
  G4double teta0 = 13.6*MeV*z*std::sqrt(t)*(1.+0.038*std::log(t))/bpc;
  return teta0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
