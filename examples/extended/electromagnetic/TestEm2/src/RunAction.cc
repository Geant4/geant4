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
/// \file electromagnetic/TestEm2/src/RunAction.cc
/// \brief Implementation of the RunAction class
//
// $Id$
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "EmAcceptance.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "G4SystemOfUnits.hh"
#include <iomanip>

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:fDet(det),fKin(kin)
{
  fRunMessenger = new RunActionMessenger(this);
  
  f_nLbin = f_nRbin = MaxBin;
    
  f_dEdL.resize(f_nLbin, 0.0);
  fSumELongit.resize(f_nLbin, 0.0);
  fSumELongitCumul.resize(f_nLbin, 0.0);
  fSumE2Longit.resize(f_nLbin, 0.0);
  fSumE2LongitCumul.resize(f_nLbin, 0.0);

  f_dEdR.resize(f_nRbin, 0.0);
  fSumERadial.resize(f_nRbin, 0.0);
  fSumERadialCumul.resize(f_nRbin, 0.0);
  fSumE2Radial.resize(f_nRbin, 0.0);
  fSumE2RadialCumul.resize(f_nRbin, 0.0);
    
  fEdeptrue = 1.;
  fRmstrue  = 1.;
  fLimittrue = DBL_MAX;

  fChargedStep = 0.0;
  fNeutralStep = 0.0;    
  
  fVerbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete fRunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BookHisto()
{
  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  // Open an output file
  //
  fHistoName[0] = "testem2";
  analysisManager->OpenFile(fHistoName[0]);    
  analysisManager->SetVerboseLevel(1);
  G4String extension = analysisManager->GetFileType();
  fHistoName[1] = fHistoName[0] + "." + extension;  

  // Creating histograms
  //
  G4double Ekin = fKin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = fDet->GetdLradl();
  G4double dRradl = fDet->GetdRradl();
  
  analysisManager->SetFirstHistoId(1);   
  analysisManager->CreateH1( "1","total energy deposit(percent of Einc)",
                                  110,0.,110.);

  analysisManager->CreateH1( "2","total charged tracklength (radl)",
                                  110,0.,110.*Ekin/GeV);

  analysisManager->CreateH1( "3","total neutral tracklength (radl)",
                                  110,0.,1100.*Ekin/GeV);

  analysisManager->CreateH1( "4","longit energy profile (% of E inc)",
                                    f_nLbin,0.,f_nLbin*dLradl);
                                    
  analysisManager->CreateH1( "5","rms on longit Edep (% of E inc)",
                                  f_nLbin,0.,f_nLbin*dLradl);

  G4double Zmin=0.5*dLradl, Zmax=Zmin+f_nLbin*dLradl;
  analysisManager->CreateH1( "6","cumul longit energy dep (% of E inc)",
                                  f_nLbin,Zmin,Zmax);
                                    
  analysisManager->CreateH1( "7","rms on cumul longit Edep (% of E inc)",
                                  f_nLbin,Zmin,Zmax);

  analysisManager->CreateH1( "8","radial energy profile (% of E inc)",
                                  f_nRbin,0.,f_nRbin*dRradl);
                                                                        
  analysisManager->CreateH1( "9","rms on radial Edep (% of E inc)",
                                  f_nRbin,0.,f_nRbin*dRradl);            

  G4double Rmin=0.5*dRradl, Rmax=Rmin+f_nRbin*dRradl;
  analysisManager->CreateH1("10","cumul radial energy dep (% of E inc)",
                                  f_nRbin,Rmin,Rmax);

  analysisManager->CreateH1("11","rms on cumul radial Edep (% of E inc)",
                                  f_nRbin,Rmin,Rmax);                    
                                    
 G4cout << "\n----> Histogram file is opened in " << fHistoName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::CleanHisto()
{
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  
  analysisManager->Write();
  analysisManager->CloseFile();

  delete G4AnalysisManager::Instance();    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  ////G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  fChargedStep = 0.0;
  fNeutralStep = 0.0;    

  f_nLbin = fDet->GetnLtot();
  f_nRbin = fDet->GetnRtot();

  //initialize arrays of cumulative energy deposition
  //
  for (G4int i=0; i<f_nLbin; i++) 
     fSumELongit[i]=fSumE2Longit[i]=fSumELongitCumul[i]=fSumE2LongitCumul[i]=0.;
  
  for (G4int j=0; j<f_nRbin; j++)
     fSumERadial[j]=fSumE2Radial[j]=fSumERadialCumul[j]=fSumE2RadialCumul[j]=0.;

  //initialize track length
  fSumChargTrLength=fSum2ChargTrLength=fSumNeutrTrLength=fSum2NeutrTrLength=0.;

  //histograms
  //
  if (aRun->GetRunID() == 0) BookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::InitializePerEvent()
{
  //initialize arrays of energy deposit per bin
  for (G4int i=0; i<f_nLbin; i++)
     { f_dEdL[i] = 0.; }
     
  for (G4int j=0; j<f_nRbin; j++)
     { f_dEdR[j] = 0.; }     
  
  //initialize tracklength 
    fChargTrLength = fNeutrTrLength = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillPerEvent()
{
  //accumulate statistic
  //
  G4double dLCumul = 0.;
  for (G4int i=0; i<f_nLbin; i++)
     {
      fSumELongit[i]  += f_dEdL[i];
      fSumE2Longit[i] += f_dEdL[i]*f_dEdL[i];
      dLCumul        += f_dEdL[i];
      fSumELongitCumul[i]  += dLCumul;
      fSumE2LongitCumul[i] += dLCumul*dLCumul;
     }

  G4double dRCumul = 0.;
  for (G4int j=0; j<f_nRbin; j++)
     {
      fSumERadial[j]  += f_dEdR[j];
      fSumE2Radial[j] += f_dEdR[j]*f_dEdR[j];
      dRCumul        += f_dEdR[j];
      fSumERadialCumul[j]  += dRCumul;
      fSumE2RadialCumul[j] += dRCumul*dRCumul;
     }

  fSumChargTrLength  += fChargTrLength;
  fSum2ChargTrLength += fChargTrLength*fChargTrLength;
  fSumNeutrTrLength  += fNeutrTrLength;
  fSum2NeutrTrLength += fNeutrTrLength*fNeutrTrLength;

  //fill histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  

  G4double Ekin=fKin->GetParticleGun()->GetParticleEnergy();
  G4double mass=fKin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=fDet->GetMaterial()->GetRadlen();

  analysisManager->FillH1(1, 100.*dLCumul/(Ekin+mass));
  analysisManager->FillH1(2, fChargTrLength/radl);
  analysisManager->FillH1(3, fNeutrTrLength/radl);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double kinEnergy = fKin->GetParticleGun()->GetParticleEnergy();
  assert(NbOfEvents*kinEnergy > 0);

  fChargedStep /= G4double(NbOfEvents);
  fNeutralStep /= G4double(NbOfEvents);    

  G4double mass=fKin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double norme = 100./(NbOfEvents*(kinEnergy+mass));
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      
  //longitudinal
  //
  G4double dLradl = fDet->GetdLradl();

  MyVector MeanELongit(f_nLbin),      rmsELongit(f_nLbin);
  MyVector MeanELongitCumul(f_nLbin), rmsELongitCumul(f_nLbin);

  G4int i;
  for (i=0; i<f_nLbin; i++)
   {
    MeanELongit[i] = norme*fSumELongit[i];
     rmsELongit[i] = norme*std::sqrt(std::fabs(NbOfEvents*fSumE2Longit[i]
                                - fSumELongit[i]*fSumELongit[i]));

    MeanELongitCumul[i] = norme*fSumELongitCumul[i];
     rmsELongitCumul[i] = norme*std::sqrt(std::fabs(NbOfEvents*fSumE2LongitCumul[i]
                                    - fSumELongitCumul[i]*fSumELongitCumul[i]));

    G4double bin = (i+0.5)*dLradl;
    analysisManager->FillH1(4, bin,MeanELongit[i]/dLradl);
    analysisManager->FillH1(5, bin, rmsELongit[i]/dLradl);      
    bin = (i+1)*dLradl;
    analysisManager->FillH1(6, bin,MeanELongitCumul[i]);
    analysisManager->FillH1(7, bin, rmsELongitCumul[i]);
   }

  //radial
  //
  G4double dRradl = fDet->GetdRradl();

  MyVector MeanERadial(f_nRbin),      rmsERadial(f_nRbin);
  MyVector MeanERadialCumul(f_nRbin), rmsERadialCumul(f_nRbin);

  for (i=0; i<f_nRbin; i++)
   {
    MeanERadial[i] = norme*fSumERadial[i];
     rmsERadial[i] = norme*std::sqrt(std::fabs(NbOfEvents*fSumE2Radial[i]
                                - fSumERadial[i]*fSumERadial[i]));

    MeanERadialCumul[i] = norme*fSumERadialCumul[i];
     rmsERadialCumul[i] = norme*std::sqrt(std::fabs(NbOfEvents*fSumE2RadialCumul[i]
                                    - fSumERadialCumul[i]*fSumERadialCumul[i]));


    G4double bin = (i+0.5)*dRradl;
    analysisManager->FillH1(8, bin,MeanERadial[i]/dRradl);
    analysisManager->FillH1(9, bin, rmsERadial[i]/dRradl);      
    bin = (i+1)*dRradl;
    analysisManager->FillH1(10, bin,MeanERadialCumul[i]);
    analysisManager->FillH1(11, bin, rmsERadialCumul[i]);
   }

  //find Moliere confinement
  //
  const G4double EMoliere = 90.;
  G4double iMoliere = 0.;
  if ((MeanERadialCumul[0]       <= EMoliere) &&
      (MeanERadialCumul[f_nRbin-1] >= EMoliere)) {
    G4int imin = 0;
    while( (imin < f_nRbin-1) && (MeanERadialCumul[imin] < EMoliere) ) imin++;
    G4double ratio = (EMoliere - MeanERadialCumul[imin]) /
                     (MeanERadialCumul[imin+1] - MeanERadialCumul[imin]);
    iMoliere = 1. + imin + ratio;
  }                       
      
  //track length
  //
  norme = 1./(NbOfEvents*(fDet->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*fSumChargTrLength;
  G4double  rmsChargTrLength = 
            norme*std::sqrt(std::fabs(NbOfEvents*fSum2ChargTrLength
                                         - fSumChargTrLength*fSumChargTrLength));

  G4double MeanNeutrTrLength = norme*fSumNeutrTrLength;
  G4double  rmsNeutrTrLength = 
            norme*std::sqrt(std::fabs(NbOfEvents*fSum2NeutrTrLength
                                         - fSumNeutrTrLength*fSumNeutrTrLength));

  //print
  //

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(2);
  
  if (fVerbose) {

    G4cout << "                 LOGITUDINAL PROFILE                   "
           << "      CUMULATIVE LOGITUDINAL PROFILE" << G4endl << G4endl;

    G4cout << "        bin   " << "           Mean         rms         "
           << "        bin "   << "           Mean      rms \n" << G4endl;

    for (i=0; i<f_nLbin; i++)
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

    for (i=0; i<f_nRbin; i++)
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
  }
  
  G4cout << "\n SUMMARY \n" << G4endl;

  G4cout << "\n"
	 << " Mean number of charged steps:  " << fChargedStep << G4endl;
  G4cout << " Mean number of neutral steps:  " << fNeutralStep 
	 << "\n" << G4endl;

  G4cout << " energy deposit : "
         << std::setw(7)  << MeanELongitCumul[f_nLbin-1] << " % E0 +- "
         << std::setw(7)  <<  rmsELongitCumul[f_nLbin-1] << " % E0" << G4endl;
  G4cout << " charged traklen: "
         << std::setw(7)  << MeanChargTrLength << " radl +- "
         << std::setw(7)  <<  rmsChargTrLength << " radl" << G4endl;
  G4cout << " neutral traklen: "
         << std::setw(7)  << MeanNeutrTrLength << " radl +- "
         << std::setw(7)  <<  rmsNeutrTrLength << " radl" << G4endl;
         
  if (iMoliere > 0. ) {
    G4double RMoliere1 = iMoliere*fDet->GetdRradl();
    G4double RMoliere2 = iMoliere*fDet->GetdRlength();          
    G4cout << "\n " << EMoliere << " % confinement: radius = "
           << RMoliere1 << " radl  ("
           << G4BestUnit( RMoliere2, "Length") << ")" << "\n" << G4endl;
  }           

  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // save histos and close AnalysisFactory
  CleanHisto();
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();

  // Acceptance
  if(fLimittrue < DBL_MAX) {
    EmAcceptance acc;
    acc.BeginOfAcceptance("Total Energy in Absorber",NbOfEvents);
    G4double e = MeanELongitCumul[f_nLbin-1]/100.;
    G4double r = rmsELongitCumul[f_nLbin-1]/100.;
    acc.EmAcceptanceGauss("Edep",NbOfEvents,e,fEdeptrue,fRmstrue,fLimittrue);
    acc.EmAcceptanceGauss("Erms",NbOfEvents,r,fRmstrue,fRmstrue,2.0*fLimittrue);
    acc.EndOfAcceptance();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetEdepAndRMS(G4ThreeVector Value)
{
  fEdeptrue = Value(0);
  fRmstrue  = Value(1);
  fLimittrue= Value(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
