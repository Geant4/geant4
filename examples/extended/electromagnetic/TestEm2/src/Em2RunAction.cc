// This code implementation is the intellectual property of
// the GEANT4 collaboration.
//
// By copying, distributing or modifying the Program (or any work
// based on the Program) you indicate your acceptance of this statement,
// and all its terms.
//
// $Id: Em2RunAction.cc,v 1.4 2000-01-21 10:56:15 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include "Em2RunAction.hh"

#include "Em2DetectorConstruction.hh"
#include "Em2PrimaryGeneratorAction.hh"
#include "Em2RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "g4std/iomanip"
#include <assert.h>

#include "Randomize.hh"

#ifndef G4NOHIST
 #include "CLHEP/Hist/HBookFile.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2RunAction::Em2RunAction(Em2DetectorConstruction*   det,
                           Em2PrimaryGeneratorAction* kin)
:Em2Det(det),Em2Kin(kin)			   
{ 
  nLbin = Em2Det->GetnLtot();
  dEdL             = *new MyVector(nLbin);
  sumELongit       = *new MyVector(nLbin);
  sumELongitCumul  = *new MyVector(nLbin);
  sumE2Longit      = *new MyVector(nLbin);
  sumE2LongitCumul = *new MyVector(nLbin);
  
  gammaFlux        = *new MyVector(nLbin);
  electronFlux     = *new MyVector(nLbin);
  positronFlux     = *new MyVector(nLbin);
    
  nRbin = Em2Det->GetnRtot();
  dEdR             = *new MyVector(nRbin);
  sumERadial       = *new MyVector(nRbin);
  sumERadialCumul  = *new MyVector(nRbin);
  sumE2Radial      = *new MyVector(nRbin);
  sumE2RadialCumul = *new MyVector(nRbin);
   
  runMessenger = new Em2RunActionMessenger(this);   
  saveRndm = 1;
  
#ifndef G4NOHIST
   hbookManager = NULL;
#endif
    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em2RunAction::~Em2RunAction()
{
  cleanHisto();    
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2RunAction::bookHisto()
{
#ifndef G4NOHIST
  hbookManager = new HBookFile("TestEm2.histo", 68);
  assert (hbookManager != 0);
  
  G4double Ekin = Em2Kin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = Em2Det->GetdLradl();
  G4double dRradl = Em2Det->GetdRradl();
  
  histo1 = hbookManager->histogram("total energy deposit (percent of E inc)",
                                    100,0.,100.);
                                    
  histo2 = hbookManager->histogram("total charged tracklength (radl)",
                                    100,0.,100.*Ekin/GeV);
                                    
  histo3 = hbookManager->histogram("total neutral tracklength (radl)",
                                    100,0.,1000.*Ekin/GeV);
                                    
  histo4 = hbookManager->histogram("longit energy profile (% of E inc)",
                                    nLbin,0.,nLbin*dLradl);
                                    
  G4double Zmin=0.5*dLradl, Zmax=Zmin+nLbin*dLradl;                                 
  histo5 = hbookManager->histogram("cumul longit energy dep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo6 = hbookManager->histogram("rms on cumul longit Edep (% of E inc)",
                                    nLbin,Zmin,Zmax);
                                    
  histo7 = hbookManager->histogram("nb of gamma per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo8 = hbookManager->histogram("nb of positron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  histo9 = hbookManager->histogram("nb of electron per plane",
                                    nLbin,Zmin,Zmax);
                                    
  hist10 = hbookManager->histogram("radial energy profile (% of E inc)",
                                    nRbin,0.,nRbin*dRradl);
                                    
  G4double Rmin=0.5*dRradl, Rmax=Rmin+nRbin*dRradl;                                 
  hist11 = hbookManager->histogram("cumul radial energy dep (% of E inc)",
                                    nRbin,Rmin,Rmax);
                                    
  hist12 = hbookManager->histogram("rms on cumul radial Edep (% of E inc)",
                                    nRbin,Rmin,Rmax);
#endif				                                                                                                                                                                                                                                                                                                                                     
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2RunAction::cleanHisto()
{
#ifndef G4NOHIST
  delete histo1; delete histo2; delete histo3;
  delete histo4; delete histo5; delete histo6;
  delete histo7; delete histo8; delete histo9;
  delete hist10; delete hist11; delete hist12;
  delete hbookManager;
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  if (saveRndm > 0)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("beginOfRun.rndm");
    }    
  
  //reshape arrays if needed  
  //
  G4bool rebin(false);
  
  G4int nLtot = Em2Det->GetnLtot();
  if(nLbin != nLtot)
    {dEdL.reshape(nLbin=nLtot);
     sumELongit.reshape(nLbin)     ; sumE2Longit.reshape(nLbin);
     sumELongitCumul.reshape(nLbin); sumE2LongitCumul.reshape(nLbin);
     gammaFlux.reshape(nLbin) ;
     electronFlux.reshape(nLbin)   ; positronFlux.reshape(nLbin);
     rebin=true;
    }
    
  G4int nRtot = Em2Det->GetnRtot();
  if(nRbin != nRtot)
    {dEdR.reshape(nRbin=nRtot);
     sumERadial.reshape(nRbin)     ; sumE2Radial.reshape(nRbin);
     sumERadialCumul.reshape(nRbin); sumE2RadialCumul.reshape(nRbin);
     rebin=true;
    }
          
  //initialize arrays of cumulative energy deposition
  //    
  for (G4int i=0; i<nLbin; i++)
     sumELongit(i)=sumE2Longit(i)=sumELongitCumul(i)=sumE2LongitCumul(i)=0.;
     
  for (G4int j=0; j<nRbin; j++)
     sumELongit(j)=sumE2Longit(j)=sumELongitCumul(j)=sumE2LongitCumul(j)=0.;
             
  //initialize track length
  sumChargTrLength=sum2ChargTrLength=sumNeutrTrLength=sum2NeutrTrLength=0.;

  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
  else if (rebin) {cleanHisto(); bookHisto();}
  
  //drawing
  // 
  if (G4VVisManager::GetConcreteInstance())
    {
      G4UImanager* UI = G4UImanager::GetUIpointer(); 
      UI->ApplyCommand("/vis/scene/notifyHandlers");
    } 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2RunAction::fillPerEvent()
{
  //accumulate statistic
  //
  G4double dLCumul = 0.;     
  for (G4int i=0; i<nLbin; i++) 
     {
      sumELongit(i)  += dEdL(i);
      sumE2Longit(i) += dEdL(i)*dEdL(i);
      dLCumul        += dEdL(i);
      sumELongitCumul(i)  += dLCumul;
      sumE2LongitCumul(i) += dLCumul*dLCumul;
     }
     
  G4double dRCumul = 0.;     
  for (G4int j=0; j<nRbin; j++) 
     {
      sumERadial(j)  += dEdR(j);
      sumE2Radial(j) += dEdR(j)*dEdR(j);
      dRCumul        += dEdR(j);
      sumERadialCumul(j)  += dRCumul;
      sumE2RadialCumul(j) += dRCumul*dRCumul;
     }
          
  sumChargTrLength  += ChargTrLength;
  sum2ChargTrLength += ChargTrLength*ChargTrLength;
  sumNeutrTrLength  += NeutrTrLength;
  sum2NeutrTrLength += NeutrTrLength*NeutrTrLength;
  
#ifndef G4NOHIST  
  //fill histograms
  //
  G4double Ekin=Em2Kin->GetParticleGun()->GetParticleEnergy();
  G4double mass=Em2Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=Em2Det->GetMaterial()->GetRadlen();
  
  histo1->accumulate(100.*dLCumul/(Ekin+mass));
  histo2->accumulate(ChargTrLength/radl);
  histo3->accumulate(NeutrTrLength/radl);
#endif      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em2RunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
     
  //compute and print statistic
  //     
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double kinEnergy = Em2Kin->GetParticleGun()->GetParticleEnergy();   
  assert(NbOfEvents*kinEnergy > 0);
  
  G4double mass=Em2Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();  
  G4double norme = 100./(NbOfEvents*(kinEnergy+mass));
  
  //longitudinal
  //
  G4double dLradl = Em2Det->GetdLradl();    
  
  MyVector MeanELongit     (nLbin), rmsELongit     (nLbin);
  MyVector MeanELongitCumul(nLbin), rmsELongitCumul(nLbin);
   
  G4int i; G4double bin;  
  for (i=0; i<nLbin; i++)
   {
    MeanELongit(i) = norme*sumELongit(i);
     rmsELongit(i) = norme*sqrt(abs(NbOfEvents*sumE2Longit(i)
                                - sumELongit(i)*sumELongit(i)));
          
    MeanELongitCumul(i) = norme*sumELongitCumul(i);
     rmsELongitCumul(i) = norme*sqrt(abs(NbOfEvents*sumE2LongitCumul(i) 
                                    - sumELongitCumul(i)*sumELongitCumul(i)));

    gammaFlux   (i) /= NbOfEvents;
    electronFlux(i) /= NbOfEvents;
    positronFlux(i) /= NbOfEvents;                                    

#ifndef G4NOHIST                                    
    bin = i*dLradl;                                
    histo4->accumulate(bin,MeanELongit(i)/dLradl);
    bin = (i+1)*dLradl;
    histo5->accumulate(bin,MeanELongitCumul(i));
    histo6->accumulate(bin, rmsELongitCumul(i));
    
    histo7->accumulate(bin, gammaFlux(i));
    histo8->accumulate(bin, electronFlux(i));
    histo9->accumulate(bin, positronFlux(i));
#endif                                            
   }
   
  //radial
  // 
  G4double dRradl = Em2Det->GetdRradl();    
  
  MyVector MeanERadial     (nRbin), rmsERadial     (nRbin);
  MyVector MeanERadialCumul(nRbin), rmsERadialCumul(nRbin);
    
  for (i=0; i<nRbin; i++)
   {
    MeanERadial(i) = norme*sumERadial(i);
     rmsERadial(i) = norme*sqrt(abs(NbOfEvents*sumE2Radial(i)
                                - sumERadial(i)*sumERadial(i)));
          
    MeanERadialCumul(i) = norme*sumERadialCumul(i);
     rmsERadialCumul(i) = norme*sqrt(abs(NbOfEvents*sumE2RadialCumul(i) 
                                    - sumERadialCumul(i)*sumERadialCumul(i)));
                                 
#ifndef G4NOHIST                                     
    bin = i*dRradl;                                
    hist10->accumulate(bin,MeanERadial(i)/dRradl);
    bin = (i+1)*dRradl;
    hist11->accumulate(bin,MeanERadialCumul(i));
    hist12->accumulate(bin, rmsERadialCumul(i));
#endif                                           
   }

  //track length
  //
  norme = 1./(NbOfEvents*(Em2Det->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*sumChargTrLength;
  G4double  rmsChargTrLength = norme*sqrt(abs(NbOfEvents*sum2ChargTrLength
                                         - sumChargTrLength*sumChargTrLength));
                                         
  G4double MeanNeutrTrLength = norme*sumNeutrTrLength;
  G4double  rmsNeutrTrLength = norme*sqrt(abs(NbOfEvents*sum2NeutrTrLength
                                         - sumNeutrTrLength*sumNeutrTrLength));                                         
   
  //print
  // 

  G4long oldform = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
  G4int  oldprec = G4cout.precision(2);
  
  G4cout << "                 LATERAL PROFILE                   "
         << "      CUMULATIVE LATERAL PROFILE" << G4endl << G4endl;  
     
  G4cout << "        bin   " << "           Mean         rms         " 
         << "        bin "   << "           Mean      rms \n" << G4endl;
                                         
  for (i=0; i<nLbin; i++)
   {
     G4double inf=i*dLradl, sup=inf+dLradl;
       
     G4cout << G4std::setw(8) << inf << "->" << G4std::setw(5) << sup << " radl: " 
                                      << G4std::setw(7) << MeanELongit(i) << "%  " 
                                      << G4std::setw(9) << rmsELongit(i) << "%       "                                  
                       << "      0->" << G4std::setw(5) << sup << " radl: " 
                                      << G4std::setw(7) << MeanELongitCumul(i) << "%  " 
                                      << G4std::setw(7) << rmsELongitCumul(i) << "% " 
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
       
     G4cout << G4std::setw(8) << inf << "->" << G4std::setw(5) << sup << " radl: " 
                                      << G4std::setw(7) << MeanERadial(i) << "%  " 
                                      << G4std::setw(9) << rmsERadial(i) << "%       "                                  
                       << "      0->" << G4std::setw(5) << sup << " radl: " 
                                      << G4std::setw(7) << MeanERadialCumul(i) << "%  " 
                                      << G4std::setw(7) << rmsERadialCumul(i) << "% " 
            <<G4endl;
   }  
  G4cout << G4endl;
  G4cout << G4std::setw(37) << "SUMMARY" << G4endl;
  G4cout << G4std::setw(42) << "energy deposit : " 
         << G4std::setw(7)  << MeanELongitCumul(nLbin-1) << " % E0 +- "
         << G4std::setw(7)  <<  rmsELongitCumul(nLbin-1) << " % E0" << G4endl;
  G4cout << G4std::setw(42) << "charged traklen: " 
         << G4std::setw(7)  << MeanChargTrLength << " radl +- "
         << G4std::setw(7)  <<  rmsChargTrLength << " radl" << G4endl;
  G4cout << G4std::setw(42) << "neutral traklen: " 
         << G4std::setw(7)  << MeanNeutrTrLength << " radl +- "
         << G4std::setw(7)  <<  rmsNeutrTrLength << " radl" << G4endl;
                 
                
  G4cout.setf(oldform,G4std::ios::floatfield);
  G4cout.precision(oldprec);

  // save Rndm status
  if (saveRndm == 1)
    { HepRandom::showEngineStatus();
      HepRandom::saveEngineStatus("endOfRun.rndm");
    }
    
#ifndef G4NOHIST   
  // Write histogram file 
  hbookManager->write();
#endif                              
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
