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
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// 08.03.01 Hisaya: Adapted MyVector for STL

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "EmAcceptance.hh"
#include "G4DataVector.hh"
#include "G4Threading.hh"

#include <iomanip>

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction*   det,
		     PrimaryGeneratorAction* kin)
  :Det(det),Kin(kin),NbOfEvents(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::Reset()
{
  NbOfEvents = 0;
  nLbin = Det->GetnLtot();
  dEdL.resize(nLbin);
  sumELongit.resize(nLbin);
  sumELongitCumul.resize(nLbin);
  sumE2Longit.resize(nLbin);
  sumE2LongitCumul.resize(nLbin);

  gammaFlux.resize(nLbin);
  electronFlux.resize(nLbin);
  positronFlux.resize(nLbin);

  nRbin = Det->GetnRtot();
  dEdR.resize(nRbin);
  sumERadial.resize(nRbin);
  sumERadialCumul.resize(nRbin);
  sumE2Radial.resize(nRbin);
  sumE2RadialCumul.resize(nRbin);

  //initialize arrays of cumulative energy deposition
  //
  for (G4int i=0; i<nLbin; i++) {
    sumELongit[i]=sumE2Longit[i]=sumELongitCumul[i]=sumE2LongitCumul[i]=0.;
    gammaFlux[i]=electronFlux[i]=positronFlux[i]=0.;
  }

  for (G4int j=0; j<nRbin; j++) {
     sumERadial[j]=sumE2Radial[j]=sumERadialCumul[j]=sumE2RadialCumul[j]=0.;
  }
  //initialize track length
  sumChargTrLength=sum2ChargTrLength=sumNeutrTrLength=sum2NeutrTrLength=0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  Reset();
  G4cout << "### Run " << aRun->GetRunID() << " start. nLbin= " << nLbin 
	 << G4endl;

  // save Rndm status
  // G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent()
{
  //accumulate statistic
  //
  G4double dLCumul = 0.;
  for (G4int i=0; i<nLbin; i++)
     {
      sumELongit[i]  += dEdL[i];
      sumE2Longit[i] += dEdL[i]*dEdL[i];
      dLCumul        += dEdL[i];
      sumELongitCumul[i]  += dLCumul;
      sumE2LongitCumul[i] += dLCumul*dLCumul;
     }

  G4double dRCumul = 0.;
  for (G4int j=0; j<nRbin; j++)
     {
      sumERadial[j]  += dEdR[j];
      sumE2Radial[j] += dEdR[j]*dEdR[j];
      dRCumul        += dEdR[j];
      sumERadialCumul[j]  += dRCumul;
      sumE2RadialCumul[j] += dRCumul*dRCumul;
     }

  sumChargTrLength  += ChargTrLength;
  sum2ChargTrLength += ChargTrLength*ChargTrLength;
  sumNeutrTrLength  += NeutrTrLength;
  sum2NeutrTrLength += NeutrTrLength*NeutrTrLength;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  //compute and print statistic
  //
  NbOfEvents = aRun->GetNumberOfEvent();

#ifdef G4MULTITHREADED
  if(G4Threading::IsWorkerThread() == true) { return; } 

  // in real MT run there is non-zero number of RunActions
  size_t nrun = runActions.size();
  if(0 < nrun) {
    NbOfEvents = 0;
    for(size_t it = 0; it < runActions.size(); ++it) { 
      RunAction* anAction = runActions[it];
      NbOfEvents += anAction->NbOfEvents;

      for (G4int i=0; i<nLbin; i++) {
	sumELongit[i]  += anAction->sumELongit[i];
	sumE2Longit[i] += anAction->sumE2Longit[i];
	sumELongitCumul[i]  += anAction->sumELongitCumul[i];
	sumE2LongitCumul[i] += anAction->sumE2LongitCumul[i];
      }

      for (G4int j=0; j<nRbin; j++) {
	sumERadial[j]  += anAction->sumERadial[j];
	sumE2Radial[j] += anAction->sumE2Radial[j];
	sumERadialCumul[j]  += anAction->sumERadialCumul[j];
	sumE2RadialCumul[j] += anAction->sumE2RadialCumul[j];
      }

      sumChargTrLength  += anAction->sumChargTrLength;
      sum2ChargTrLength += anAction->sum2ChargTrLength;
      sumNeutrTrLength  += anAction->sumNeutrTrLength;
      sum2NeutrTrLength += anAction->sum2NeutrTrLength;
    }
  }
#endif

  G4double kinEnergy = Kin->GetParticleGun()->GetParticleEnergy();
  assert(NbOfEvents*kinEnergy > 0);

  G4double mass=Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double norme = 100./(NbOfEvents*(kinEnergy+mass));

  //longitudinal
  //
  G4double dLradl = Det->GetdLradl();

  G4DataVector MeanELongit, rmsELongit, MeanELongitCumul, rmsELongitCumul;
  MeanELongit.resize(nLbin);
  rmsELongit.resize(nLbin);
  MeanELongitCumul.resize(nLbin);
  rmsELongitCumul.resize(nLbin);
      
  G4int i;
  for (i=0; i<nLbin; i++)
    {
      MeanELongit[i] = norme*sumELongit[i];
      rmsELongit[i] = norme*std::sqrt(std::abs(NbOfEvents*sumE2Longit[i]
						- sumELongit[i]*sumELongit[i]));
      
      MeanELongitCumul[i] = norme*sumELongitCumul[i];
      rmsELongitCumul[i] = norme*std::sqrt(std::abs(NbOfEvents*sumE2LongitCumul[i]
					   - sumELongitCumul[i]*sumELongitCumul[i]));
      
      gammaFlux   [i] /= NbOfEvents;
      electronFlux[i] /= NbOfEvents;
      positronFlux[i] /= NbOfEvents;
    }
  
  //radial
  //
  G4double dRradl = Det->GetdRradl();
  
  G4DataVector MeanERadial, rmsERadial, MeanERadialCumul, rmsERadialCumul;
  MeanERadial.resize(nRbin);
  rmsERadial.resize(nRbin);
  MeanERadialCumul.resize(nRbin);
  rmsERadialCumul.resize(nRbin);
  
  for (i=0; i<nRbin; i++)
    {
      MeanERadial[i] = norme*sumERadial[i];
      rmsERadial[i] = norme*std::sqrt(std::abs(NbOfEvents*sumE2Radial[i]
						- sumERadial[i]*sumERadial[i]));
      
      MeanERadialCumul[i] = norme*sumERadialCumul[i];
      rmsERadialCumul[i] = norme*std::sqrt(std::abs(NbOfEvents*sumE2RadialCumul[i]
						     - sumERadialCumul[i]*sumERadialCumul[i]));
      
    }
  
  //track length
  //
  norme = 1./(NbOfEvents*(Det->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*sumChargTrLength;
  G4double  rmsChargTrLength = 
    norme*std::sqrt(std::abs(NbOfEvents*sum2ChargTrLength
			     - sumChargTrLength*sumChargTrLength));
  
  G4double MeanNeutrTrLength = norme*sumNeutrTrLength;
  G4double  rmsNeutrTrLength = 
    norme*std::sqrt(std::abs(NbOfEvents*sum2NeutrTrLength
			     - sumNeutrTrLength*sumNeutrTrLength));
  
  //print
  //

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(2);
  
  G4cout << "                 LATERAL PROFILE                   "
	 << "      CUMULATIVE LATERAL PROFILE" << G4endl << G4endl;
  
  G4cout << "        bin   " << "           Mean         rms         "
	 << "        bin "   << "           Mean      rms \n" << G4endl;
  
  for (i=0; i<nLbin; i++)
    {
      G4double inf=i*dLradl, sup=inf+dLradl;
      
      G4cout << std::setw(8) << inf << "->"
	     << std::setw(5) << sup << " radl: "
	     << std::setw(7) << MeanELongit[i] << "%  "
	     << std::setw(9) << rmsELongit[i] << "%       "
	     << "      0->" << std::setw(5) << sup << " radl: "
	     << std::setw(7) << MeanELongitCumul[i] << "%  "
	     << std::setw(7) << rmsELongitCumul[i] << "% "
	     <<G4endl;
    }
  
  G4cout << G4endl << G4endl << G4endl;
  
  G4cout << "                  RADIAL PROFILE                   "
	 << "      CUMULATIVE  RADIAL PROFILE" << G4endl << G4endl;
  
  G4cout << "        bin   " << "           Mean         rms         "
	 << "        bin "   << "           Mean      rms \n" << G4endl;
  
  for (i=0; i<nRbin; i++)
    {
      G4double inf=i*dRradl, sup=inf+dRradl;
      
      G4cout << std::setw(8) << inf << "->"
	     << std::setw(5) << sup << " radl: "
	     << std::setw(7) << MeanERadial[i] << "%  "
	     << std::setw(9) << rmsERadial[i] << "%       "
	     << "      0->" << std::setw(5) << sup << " radl: "
	     << std::setw(7) << MeanERadialCumul[i] << "%  "
	     << std::setw(7) << rmsERadialCumul[i] << "% "
	     <<G4endl;
    }
  G4cout << G4endl;
  G4cout << std::setw(37) << "SUMMARY" << G4endl;
  G4cout << std::setw(42) << "energy deposit : "
	 << std::setw(7)  << MeanELongitCumul[nLbin-1] << " % E0 +- "
	 << std::setw(7)  <<  rmsELongitCumul[nLbin-1] << " % E0" << G4endl;
  G4cout << std::setw(42) << "charged traklen: "
	 << std::setw(7)  << MeanChargTrLength << " radl +- "
	 << std::setw(7)  <<  rmsChargTrLength << " radl" << G4endl;
  G4cout << std::setw(42) << "neutral traklen: "
	 << std::setw(7)  << MeanNeutrTrLength << " radl +- "
	 << std::setw(7)  <<  rmsNeutrTrLength << " radl" << G4endl;
  
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
  
  // Acceptance
  G4double ltrue = Det->GetLimitEdep();
  if(ltrue < DBL_MAX) {
    EmAcceptance acc;
    acc.BeginOfAcceptance("Total Energy in Absorber",NbOfEvents);
    G4double etrue = Det->GetAverageEdep();
    G4double rtrue = Det->GetRMSEdep();
    G4double e = MeanELongitCumul[nLbin-1]/100.;
    //G4double r = rmsELongitCumul[nLbin-1]/100.;
    acc.EmAcceptanceGauss("Edep",NbOfEvents,e,etrue,rtrue,ltrue);
    //    acc.EmAcceptanceGauss("Erms",NbOfEvents,r,rtrue,rtrue,2.0*ltrue);
    acc.EndOfAcceptance();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
