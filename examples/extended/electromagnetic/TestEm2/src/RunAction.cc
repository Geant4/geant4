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
// $Id: RunAction.cc,v 1.26 2010-11-09 21:25:15 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include <iomanip>

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
:Det(det),Kin(kin)
{
  runMessenger = new RunActionMessenger(this);
  
  nLbin = nRbin = MaxBin;
    
  dEdL.resize(nLbin, 0.0);
  sumELongit.resize(nLbin, 0.0);
  sumELongitCumul.resize(nLbin, 0.0);
  sumE2Longit.resize(nLbin, 0.0);
  sumE2LongitCumul.resize(nLbin, 0.0);

  dEdR.resize(nRbin, 0.0);
  sumERadial.resize(nRbin, 0.0);
  sumERadialCumul.resize(nRbin, 0.0);
  sumE2Radial.resize(nRbin, 0.0);
  sumE2RadialCumul.resize(nRbin, 0.0);
    
  edeptrue = 1.;
  rmstrue  = 1.;
  limittrue = DBL_MAX;
  
  verbose = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
  // Create analysis manager
  // The choice of analysis technology is done via selection of a namespace
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
    
  // Open an output file
  //
  histoName[0] = "testem2";
  analysisManager->OpenFile(histoName[0]);    
  analysisManager->SetVerboseLevel(1);
  G4String extension = analysisManager->GetFileType();
  histoName[1] = histoName[0] + "." + extension;  

  // Creating histograms
  //
  G4double Ekin = Kin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = Det->GetdLradl();
  G4double dRradl = Det->GetdRradl();
  
  analysisManager->SetFirstHistoId(1);   
  analysisManager->CreateH1( "1","total energy deposit(percent of Einc)",
                                  110,0.,110.);

  analysisManager->CreateH1( "2","total charged tracklength (radl)",
                                  110,0.,110.*Ekin/GeV);

  analysisManager->CreateH1( "3","total neutral tracklength (radl)",
                                  110,0.,1100.*Ekin/GeV);

  analysisManager->CreateH1( "4","longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);
				    
  analysisManager->CreateH1( "5","rms on longit Edep (% of E inc)",
                                  nLbin,0.,nLbin*dLradl);

  G4double Zmin=0.5*dLradl, Zmax=Zmin+nLbin*dLradl;
  analysisManager->CreateH1( "6","cumul longit energy dep (% of E inc)",
                                  nLbin,Zmin,Zmax);
				    
  analysisManager->CreateH1( "7","rms on cumul longit Edep (% of E inc)",
                                  nLbin,Zmin,Zmax);

  analysisManager->CreateH1( "8","radial energy profile (% of E inc)",
                                  nRbin,0.,nRbin*dRradl);
				    				    
  analysisManager->CreateH1( "9","rms on radial Edep (% of E inc)",
                                  nRbin,0.,nRbin*dRradl);	    

  G4double Rmin=0.5*dRradl, Rmax=Rmin+nRbin*dRradl;
  analysisManager->CreateH1("10","cumul radial energy dep (% of E inc)",
                                  nRbin,Rmin,Rmax);

  analysisManager->CreateH1("11","rms on cumul radial Edep (% of E inc)",
                                  nRbin,Rmin,Rmax);		    
				    
 G4cout << "\n----> Histogram file is opened in " << histoName[1] << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::cleanHisto()
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
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();

  nLbin = Det->GetnLtot();
  nRbin = Det->GetnRtot();

  //initialize arrays of cumulative energy deposition
  //
  for (G4int i=0; i<nLbin; i++) 
     sumELongit[i]=sumE2Longit[i]=sumELongitCumul[i]=sumE2LongitCumul[i]=0.;
  
  for (G4int j=0; j<nRbin; j++)
     sumERadial[j]=sumE2Radial[j]=sumERadialCumul[j]=sumE2RadialCumul[j]=0.;

  //initialize track length
  sumChargTrLength=sum2ChargTrLength=sumNeutrTrLength=sum2NeutrTrLength=0.;

  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
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

  //fill histograms
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();  

  G4double Ekin=Kin->GetParticleGun()->GetParticleEnergy();
  G4double mass=Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=Det->GetMaterial()->GetRadlen();

  analysisManager->FillH1(1, 100.*dLCumul/(Ekin+mass));
  analysisManager->FillH1(2, ChargTrLength/radl);
  analysisManager->FillH1(3, NeutrTrLength/radl);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double kinEnergy = Kin->GetParticleGun()->GetParticleEnergy();
  assert(NbOfEvents*kinEnergy > 0);

  G4double mass=Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double norme = 100./(NbOfEvents*(kinEnergy+mass));
  
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
      
  //longitudinal
  //
  G4double dLradl = Det->GetdLradl();

  MyVector MeanELongit(nLbin),      rmsELongit(nLbin);
  MyVector MeanELongitCumul(nLbin), rmsELongitCumul(nLbin);

  G4int i;
  for (i=0; i<nLbin; i++)
   {
    MeanELongit[i] = norme*sumELongit[i];
     rmsELongit[i] = norme*std::sqrt(std::fabs(NbOfEvents*sumE2Longit[i]
                                - sumELongit[i]*sumELongit[i]));

    MeanELongitCumul[i] = norme*sumELongitCumul[i];
     rmsELongitCumul[i] = norme*std::sqrt(std::fabs(NbOfEvents*sumE2LongitCumul[i]
                                    - sumELongitCumul[i]*sumELongitCumul[i]));

    G4double bin = (i+0.5)*dLradl;
    analysisManager->FillH1(4, bin,MeanELongit[i]/dLradl);
    analysisManager->FillH1(5, bin, rmsELongit[i]/dLradl);      
    bin = (i+1)*dLradl;
    analysisManager->FillH1(6, bin,MeanELongitCumul[i]);
    analysisManager->FillH1(7, bin, rmsELongitCumul[i]);
   }

  //radial
  //
  G4double dRradl = Det->GetdRradl();

  MyVector MeanERadial(nRbin),      rmsERadial(nRbin);
  MyVector MeanERadialCumul(nRbin), rmsERadialCumul(nRbin);

  for (i=0; i<nRbin; i++)
   {
    MeanERadial[i] = norme*sumERadial[i];
     rmsERadial[i] = norme*std::sqrt(std::fabs(NbOfEvents*sumE2Radial[i]
                                - sumERadial[i]*sumERadial[i]));

    MeanERadialCumul[i] = norme*sumERadialCumul[i];
     rmsERadialCumul[i] = norme*std::sqrt(std::fabs(NbOfEvents*sumE2RadialCumul[i]
                                    - sumERadialCumul[i]*sumERadialCumul[i]));


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
      (MeanERadialCumul[nRbin-1] >= EMoliere)) {
    G4int imin = 0;
    while( (imin < nRbin-1) && (MeanERadialCumul[imin] < EMoliere) ) imin++;
    G4double ratio = (EMoliere - MeanERadialCumul[imin]) /
                     (MeanERadialCumul[imin+1] - MeanERadialCumul[imin]);
    iMoliere = 1. + imin + ratio;
  }  		     
      
  //track length
  //
  norme = 1./(NbOfEvents*(Det->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*sumChargTrLength;
  G4double  rmsChargTrLength = 
            norme*std::sqrt(std::fabs(NbOfEvents*sum2ChargTrLength
                                         - sumChargTrLength*sumChargTrLength));

  G4double MeanNeutrTrLength = norme*sumNeutrTrLength;
  G4double  rmsNeutrTrLength = 
            norme*std::sqrt(std::fabs(NbOfEvents*sum2NeutrTrLength
                                         - sumNeutrTrLength*sumNeutrTrLength));

  //print
  //

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(2);
  
  if (verbose) {

    G4cout << "                 LOGITUDINAL PROFILE                   "
           << "      CUMULATIVE LOGITUDINAL PROFILE" << G4endl << G4endl;

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
  }
  
  G4cout << "\n SUMMARY \n" << G4endl;
  G4cout << " energy deposit : "
         << std::setw(7)  << MeanELongitCumul[nLbin-1] << " % E0 +- "
         << std::setw(7)  <<  rmsELongitCumul[nLbin-1] << " % E0" << G4endl;
  G4cout << " charged traklen: "
         << std::setw(7)  << MeanChargTrLength << " radl +- "
         << std::setw(7)  <<  rmsChargTrLength << " radl" << G4endl;
  G4cout << " neutral traklen: "
         << std::setw(7)  << MeanNeutrTrLength << " radl +- "
         << std::setw(7)  <<  rmsNeutrTrLength << " radl" << G4endl;
	 
  if (iMoliere > 0. ) {
    G4double RMoliere1 = iMoliere*Det->GetdRradl();
    G4double RMoliere2 = iMoliere*Det->GetdRlength(); 	 
    G4cout << "\n " << EMoliere << " % confinement: radius = "
           << RMoliere1 << " radl  ("
	   << G4BestUnit( RMoliere2, "Length") << ")" << G4endl;
  }	   

  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // save histos and close AnalysisFactory
  cleanHisto();
  
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();

  // Acceptance
  if(limittrue < DBL_MAX) {
    EmAcceptance acc;
    acc.BeginOfAcceptance("Total Energy in Absorber",NbOfEvents);
    G4double e = MeanELongitCumul[nLbin-1]/100.;
    G4double r = rmsELongitCumul[nLbin-1]/100.;
    acc.EmAcceptanceGauss("Edep",NbOfEvents,e,edeptrue,rmstrue,limittrue);
    acc.EmAcceptanceGauss("Erms",NbOfEvents,r,rmstrue,rmstrue,2.0*limittrue);
    acc.EndOfAcceptance();
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetEdepAndRMS(G4ThreeVector Value)
{
  edeptrue = Value(0);
  rmstrue  = Value(1);
  limittrue= Value(2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
