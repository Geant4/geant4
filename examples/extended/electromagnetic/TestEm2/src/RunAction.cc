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
//
// $Id: RunAction.cc,v 1.2 2004-03-15 14:26:00 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
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

#include <iomanip>
#include <assert.h>

#include "Randomize.hh"

#ifdef USE_AIDA
 #include "AIDA/AIDA.h"
#endif

#ifdef USE_ROOT
 #include "TFile.h"
 #include "TH1F.h"
#endif


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction*   det,
                           PrimaryGeneratorAction* kin)
:Det(det),Kin(kin)			   
{ 
  nLbin = Det->GetnLtot();
  dEdL.resize(nLbin, 0.0);
  sumELongit.resize(nLbin, 0.0);
  sumELongitCumul.resize(nLbin, 0.0); 
  sumE2Longit.resize(nLbin, 0.0);     
  sumE2LongitCumul.resize(nLbin, 0.0);
  
  gammaFlux.resize(nLbin, 0.0);
  electronFlux.resize(nLbin, 0.0);
  positronFlux.resize(nLbin, 0.0);
    
  nRbin = Det->GetnRtot();
  dEdR.resize(nRbin, 0.0);
  sumERadial.resize(nRbin, 0.0);
  sumERadialCumul.resize(nRbin, 0.0);
  sumE2Radial.resize(nRbin, 0.0);
  sumE2RadialCumul.resize(nRbin, 0.0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  cleanHisto();    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
#ifdef USE_AIDA
 // Creating the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 AIDA::ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create("testem2.paw", "hbook", readOnly, createNew);

 // Creating a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);  
  G4double Ekin = Kin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = Det->GetdLradl();
  G4double dRradl = Det->GetdRradl();
  
  histo[0] = hf->createHistogram1D( "1","total energy deposit(percent of Einc)",
                                    100,0.,100.);
                                    
  histo[1] = hf->createHistogram1D( "2","total charged tracklength (radl)",
                                    100,0.,100.*Ekin/GeV);
                                    
  histo[2] = hf->createHistogram1D( "3","total neutral tracklength (radl)",
                                    100,0.,1000.*Ekin/GeV);
                                    
  histo[3] = hf->createHistogram1D( "4","longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);
                                    
  G4double Zmin=0.5*dLradl, Zmax=Zmin+nLbin*dLradl;
  histo[4] = hf->createHistogram1D( "5","cumul longit energy dep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo[5] = hf->createHistogram1D( "6","rms on cumul longit Edep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo[6] = hf->createHistogram1D( "7","nb of gamma per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[7] = hf->createHistogram1D( "8","nb of positron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[8] = hf->createHistogram1D( "9","nb of electron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[9] = hf->createHistogram1D("10","radial energy profile (% of E inc)",
                                    nRbin,0.,nRbin*dRradl);
                                    
  G4double Rmin=0.5*dRradl, Rmax=Rmin+nRbin*dRradl;
  histo[10]= hf->createHistogram1D("11","cumul radial energy dep (% of E inc)",
                                    nRbin,Rmin,Rmax);
                                    
  histo[11]= hf->createHistogram1D("12","rms on cumul radial Edep (% of E inc)",
                                    nRbin,Rmin,Rmax);
				    
  delete hf;
  delete tf;
  delete af;				    
#endif

#ifdef USE_ROOT
 // Create a ROOT file
 tree = new TFile("testem2.root","recreate");
 
 // Create the histograms
  G4double Ekin = Kin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = Det->GetdLradl();
  G4double dRradl = Det->GetdRradl();
  
  histo[0] = new TH1F( "1","total energy deposit(percent of Einc)",
                                    100,0.,100.);
                                    
  histo[1] = new TH1F( "2","total charged tracklength (radl)",
                                    100,0.,100.*Ekin/GeV);
                                    
  histo[2] = new TH1F( "3","total neutral tracklength (radl)",
                                    100,0.,1000.*Ekin/GeV);
                                    
  histo[3] = new TH1F( "4","longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);
                                    
  G4double Zmin=0.5*dLradl, Zmax=Zmin+nLbin*dLradl;
  histo[4] = new TH1F( "5","cumul longit energy dep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo[5] = new TH1F( "6","rms on cumul longit Edep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo[6] = new TH1F( "7","nb of gamma per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[7] = new TH1F( "8","nb of positron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[8] = new TH1F( "9","nb of electron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo[9] = new TH1F("10","radial energy profile (% of E inc)",
                                    nRbin,0.,nRbin*dRradl);
                                    
  G4double Rmin=0.5*dRradl, Rmax=Rmin+nRbin*dRradl;
  histo[10]= new TH1F("11","cumul radial energy dep (% of E inc)",
                                    nRbin,Rmin,Rmax);
                                    
  histo[11]= new TH1F("12","rms on cumul radial Edep (% of E inc)",
                                    nRbin,Rmin,Rmax);
#endif    
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::cleanHisto()
{
#ifdef USE_AIDA
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
  delete tree;
#endif 
  
#ifdef USE_ROOT
  tree->Write();        // Writing the histograms to the file
  tree->Close();        // and closing the file  
  delete tree;
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  //reshape arrays if needed  
  //
  G4bool rebin(false);
  
  G4int nLtot = Det->GetnLtot();
  if(nLbin != nLtot)
    {dEdL.resize(nLbin=nLtot, 0.0);
     sumELongit.resize(nLbin, 0.0);
     sumE2Longit.resize(nLbin, 0.0);
     sumELongitCumul.resize(nLbin, 0.0);
     sumE2LongitCumul.resize(nLbin, 0.0);
     gammaFlux.resize(nLbin, 0.0);
     electronFlux.resize(nLbin, 0.0);
     positronFlux.resize(nLbin, 0.0);
     rebin=true;
    }
    
  G4int nRtot = Det->GetnRtot();
  if(nRbin != nRtot)
    {dEdR.resize(nRbin=nRtot, 0.0);
     sumERadial.resize(nRbin, 0.0);
     sumE2Radial.resize(nRbin, 0.0);
     sumERadialCumul.resize(nRbin, 0.0);
     sumE2RadialCumul.resize(nRbin, 0.0);
     rebin=true;
    }
          
  //initialize arrays of cumulative energy deposition
  //    
  for (G4int i=0; i<nLbin; i++)
     sumELongit[i]=sumE2Longit[i]=sumELongitCumul[i]=sumE2LongitCumul[i]=0.;
     
  for (G4int j=0; j<nRbin; j++)
     sumELongit[j]=sumE2Longit[j]=sumELongitCumul[j]=sumE2LongitCumul[j]=0.;
             
  //initialize track length
  sumChargTrLength=sum2ChargTrLength=sumNeutrTrLength=sum2NeutrTrLength=0.;

  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
  else if (rebin) {cleanHisto(); bookHisto();}
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
  
#ifdef USE_AIDA  
  //fill histograms
  //
  G4double Ekin=Kin->GetParticleGun()->GetParticleEnergy();
  G4double mass=Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=Det->GetMaterial()->GetRadlen();
  
  histo[0]->fill(100.*dLCumul/(Ekin+mass));
  histo[1]->fill(ChargTrLength/radl);
  histo[2]->fill(NeutrTrLength/radl);
#endif

#ifdef USE_ROOT  
  //fill histograms
  //
  G4double Ekin=Kin->GetParticleGun()->GetParticleEnergy();
  G4double mass=Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=Det->GetMaterial()->GetRadlen();
  
  histo[0]->Fill(100.*dLCumul/(Ekin+mass));
  histo[1]->Fill(ChargTrLength/radl);
  histo[2]->Fill(NeutrTrLength/radl);
#endif     
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
  
  //longitudinal
  //
  G4double dLradl = Det->GetdLradl();    
  
  MyVector MeanELongit(nLbin),      rmsELongit(nLbin);
  MyVector MeanELongitCumul(nLbin), rmsELongitCumul(nLbin);
   
  G4int i;   
  for (i=0; i<nLbin; i++)
   {
    MeanELongit[i] = norme*sumELongit[i];
     rmsELongit[i] = norme*sqrt(abs(NbOfEvents*sumE2Longit[i]
                                - sumELongit[i]*sumELongit[i]));
          
    MeanELongitCumul[i] = norme*sumELongitCumul[i];
     rmsELongitCumul[i] = norme*sqrt(abs(NbOfEvents*sumE2LongitCumul[i] 
                                    - sumELongitCumul[i]*sumELongitCumul[i]));

    gammaFlux   [i] /= NbOfEvents;
    electronFlux[i] /= NbOfEvents;
    positronFlux[i] /= NbOfEvents;                                    

#ifdef USE_AIDA                                    
    G4double bin = i*dLradl;                                
    histo[3]->fill(bin,MeanELongit[i]/dLradl);
    bin = (i+1)*dLradl;
    histo[4]->fill(bin,MeanELongitCumul[i]);
    histo[5]->fill(bin, rmsELongitCumul[i]);
    
    histo[6]->fill(bin, gammaFlux[i]);
    histo[7]->fill(bin, positronFlux[i]);
    histo[8]->fill(bin, electronFlux[i]);
#endif

#ifdef USE_ROOT                                    
    G4double bin = i*dLradl;                                
    histo[3]->Fill(bin,MeanELongit[i]/dLradl);
    bin = (i+1)*dLradl;
    histo[4]->Fill(bin,MeanELongitCumul[i]);
    histo[5]->Fill(bin, rmsELongitCumul[i]);
    
    histo[6]->Fill(bin, gammaFlux[i]);
    histo[7]->Fill(bin, positronFlux[i]);
    histo[8]->Fill(bin, electronFlux[i]);
#endif                                           
   }
   
  //radial
  // 
  G4double dRradl = Det->GetdRradl();    
  
  MyVector MeanERadial(nRbin),      rmsERadial(nRbin);
  MyVector MeanERadialCumul(nRbin), rmsERadialCumul(nRbin);
    
  for (i=0; i<nRbin; i++)
   {
    MeanERadial[i] = norme*sumERadial[i];
     rmsERadial[i] = norme*sqrt(abs(NbOfEvents*sumE2Radial[i]
                                - sumERadial[i]*sumERadial[i]));
          
    MeanERadialCumul[i] = norme*sumERadialCumul[i];
     rmsERadialCumul[i] = norme*sqrt(abs(NbOfEvents*sumE2RadialCumul[i] 
                                    - sumERadialCumul[i]*sumERadialCumul[i]));
                                 
#ifdef USE_AIDA                                     
    G4double bin = i*dRradl;                                
    histo[ 9]->fill(bin,MeanERadial[i]/dRradl);
    bin = (i+1)*dRradl;
    histo[10]->fill(bin,MeanERadialCumul[i]);
    histo[11]->fill(bin, rmsERadialCumul[i]);
#endif

#ifdef USE_ROOT                                     
    G4double bin = i*dRradl;                                
    histo[ 9]->Fill(bin,MeanERadial[i]/dRradl);
    bin = (i+1)*dRradl;
    histo[10]->Fill(bin,MeanERadialCumul[i]);
    histo[11]->Fill(bin, rmsERadialCumul[i]);
#endif                                           
   }

  //track length
  //
  norme = 1./(NbOfEvents*(Det->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*sumChargTrLength;
  G4double  rmsChargTrLength = norme*sqrt(abs(NbOfEvents*sum2ChargTrLength
                                         - sumChargTrLength*sumChargTrLength));
                                         
  G4double MeanNeutrTrLength = norme*sumNeutrTrLength;
  G4double  rmsNeutrTrLength = norme*sqrt(abs(NbOfEvents*sum2NeutrTrLength
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
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
