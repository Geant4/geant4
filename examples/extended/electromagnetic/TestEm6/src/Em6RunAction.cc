// Em6RunAction.cc

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em6RunAction.hh"

#include "Em6DetectorConstruction.hh"
#include "Em6PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "G4MuonPlus.hh" // for muon mass

#include "g4std/iomanip"
#include <assert.h>

#include "Randomize.hh"
// #include "cmath" // for fabs

#ifndef G4NOHIST
  #include "CLHEP/Hist/HBookFile.h"
  #include "CLHEP/Hist/TupleManager.h"
  HepTuple *heptuple1;
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6RunAction::Em6RunAction(Em6DetectorConstruction*   det,
                           Em6PrimaryGeneratorAction* kin)
:Em6Det(det),Em6Kin(kin),NGamMuMu(0)
{
  nLbin = Em6Det->GetnLtot();
  dEdL.resize(nLbin, 0.0);
  sumELongit.resize(nLbin, 0.0);
  sumELongitCumul.resize(nLbin, 0.0);
  sumE2Longit.resize(nLbin, 0.0);
  sumE2LongitCumul.resize(nLbin, 0.0);

  gammaFlux.resize(nLbin, 0.0);
  electronFlux.resize(nLbin, 0.0);
  positronFlux.resize(nLbin, 0.0);

#ifndef G4NOHIST
   hbookManager = NULL;
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em6RunAction::~Em6RunAction()
{
  cleanHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6RunAction::bookHisto()
{
#ifndef G4NOHIST
  hbookManager = new HBookFile("TestEm6.rzfile", 68);
  assert (hbookManager != 0);

  G4double Ekin = Em6Kin->GetParticleGun()->GetParticleEnergy();
  G4double dLradl = Em6Det->GetdLradl();

  histo1 = hbookManager->histogram("total energy deposit (percent of E inc)",
                                    100,0.,0.1);

  histo2 = hbookManager->histogram("total charged tracklength (radl)",
                                    100,0.,0.001*Ekin/GeV);

  histo3 = hbookManager->histogram("total neutral tracklength (radl)",
                                    100,0.,1.*Ekin/GeV);

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

  hist10 = hbookManager->histogram("1/(1+(theta+[g]+)**2)",100,0,1.);
  hist11 = hbookManager->histogram("log10(theta+ [g]+)",100,-3.,1.);
  hist12 = hbookManager->histogram("log10(theta- [g]-)",100,-3.,1.);
  hist13 = hbookManager->histogram("log10(theta+ [g]+ -theta- [g]-)",
            100,-3.,1.);
  hist14 = hbookManager->histogram("xPlus",100,0.,1.);

  heptuple1=hbookManager->ntuple("MuPair");
#endif
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6RunAction::cleanHisto()
{
#ifndef G4NOHIST
  // Write histogram file
  hbookManager->write();

  delete histo1; delete histo2; delete histo3;
  delete histo4; delete histo5; delete histo6;
  delete histo7; delete histo8; delete histo9;
  delete hist10; delete hist11; delete hist12;
  delete hist13; delete hist14;
  delete hbookManager;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  //reshape arrays if needed
  //
  G4bool rebin(false);

  G4int nLtot = Em6Det->GetnLtot();
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

  //initialize arrays of cumulative energy deposition
  //
  for (G4int i=0; i<nLbin; i++)
     sumELongit[i]=sumE2Longit[i]=sumELongitCumul[i]=sumE2LongitCumul[i]=0.;

  //initialize track length
  sumChargTrLength=sum2ChargTrLength=sumNeutrTrLength=sum2NeutrTrLength=0.;

  //histograms
  //
  if (aRun->GetRunID() == 0) bookHisto();
  else if (rebin) {cleanHisto(); bookHisto();}

  //drawing
  //
  if (G4VVisManager::GetConcreteInstance())
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6RunAction::fillPerEvent()
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

  sumChargTrLength  += ChargTrLength;
  sum2ChargTrLength += ChargTrLength*ChargTrLength;
  sumNeutrTrLength  += NeutrTrLength;
  sum2NeutrTrLength += NeutrTrLength*NeutrTrLength;

#ifndef G4NOHIST
  static G4double Mmuon=G4MuonPlus::MuonPlus()->GetPDGMass();
  //fill histograms
  //
  G4double Ekin=Em6Kin->GetParticleGun()->GetParticleEnergy();
  G4double mass=Em6Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double radl=Em6Det->GetMaterial()->GetRadlen();

  histo1->accumulate(100.*dLCumul/(Ekin+mass));
  histo2->accumulate(ChargTrLength/radl);
  histo3->accumulate(NeutrTrLength/radl);
  if(Egam>0) // there was a gamma -> mu+mu- interaction
  { heptuple1->column("Egamma",Egam/GeV);
    heptuple1->column("xPlus",xPlus);
    heptuple1->column("thetaPlus",thetaPlus);
    heptuple1->column("thetaMinus",thetaMinus);
    heptuple1->dumpData();
    
    G4double GammaPlus=Egam*xPlus/Mmuon;
    hist10->accumulate(1./(1.+pow(thetaPlus*GammaPlus,2)));
    hist11->accumulate(log10(thetaPlus*GammaPlus));

    G4double GammaMinus=Egam*(1.-xPlus)/Mmuon;
    hist12->accumulate(log10(thetaMinus*GammaMinus));

    hist13->accumulate(log10(fabs(thetaPlus*GammaPlus-thetaMinus*GammaMinus)));
    hist14->accumulate(xPlus);
  }
#endif

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em6RunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double kinEnergy = Em6Kin->GetParticleGun()->GetParticleEnergy();
  assert(NbOfEvents*kinEnergy > 0);

  G4double mass=Em6Kin->GetParticleGun()->GetParticleDefinition()->GetPDGMass();
  G4double norme = 100./(NbOfEvents*(kinEnergy+mass));

  //longitudinal
  //
  G4double dLradl = Em6Det->GetdLradl();

  MyVector MeanELongit(nLbin),      rmsELongit(nLbin);
  MyVector MeanELongitCumul(nLbin), rmsELongitCumul(nLbin);

  G4int i; G4double bin;
  for (i=0; i<nLbin; i++)
   {
    MeanELongit[i] = norme*sumELongit[i];
     rmsELongit[i] = norme*sqrt(fabs(NbOfEvents*sumE2Longit[i]
                                - sumELongit[i]*sumELongit[i]));

    MeanELongitCumul[i] = norme*sumELongitCumul[i];
     rmsELongitCumul[i] = norme*sqrt(fabs(NbOfEvents*sumE2LongitCumul[i]
                                    - sumELongitCumul[i]*sumELongitCumul[i]));

    gammaFlux   [i] /= NbOfEvents;
    electronFlux[i] /= NbOfEvents;
    positronFlux[i] /= NbOfEvents;

#ifndef G4NOHIST
    bin = i*dLradl;
    histo4->accumulate(bin,MeanELongit[i]/dLradl);
    bin = (i+1)*dLradl;
    histo5->accumulate(bin,MeanELongitCumul[i]);
    histo6->accumulate(bin, rmsELongitCumul[i]);

    histo7->accumulate(bin, gammaFlux[i]);
    histo8->accumulate(bin, positronFlux[i]);
    histo9->accumulate(bin, electronFlux[i]);
#endif
   }


  //track length
  //
  norme = 1./(NbOfEvents*(Em6Det->GetMaterial()->GetRadlen()));
  G4double MeanChargTrLength = norme*sumChargTrLength;
  G4double  rmsChargTrLength = norme*sqrt(fabs(NbOfEvents*sum2ChargTrLength
                                         - sumChargTrLength*sumChargTrLength));

  G4double MeanNeutrTrLength = norme*sumNeutrTrLength;
  G4double  rmsNeutrTrLength = norme*sqrt(fabs(NbOfEvents*sum2NeutrTrLength
                                         - sumNeutrTrLength*sumNeutrTrLength));

  //print
  //
  G4cout << "There were " << NGamMuMu << " gamma->mu+mu- conversions" << '\n';

  G4long oldform = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
  G4int  oldprec = G4cout.precision(2);

  G4cout << "                 LATERAL PROFILE                   "
         << "      CUMULATIVE LATERAL PROFILE" << G4endl << G4endl;

  G4cout << "        bin    " << "           Mean         rms         "
         << "        bin  "   << "           Mean      rms \n" << G4endl;

  for (i=0; i<nLbin; i++)
   {
     G4double inf=i*dLradl, sup=inf+dLradl;

     G4cout << G4std::setw(8) << inf << "->"
            << G4std::setw(6) << sup << " radl: "
            << setprecision(4)
            << G4std::setw(7) << MeanELongit[i] << "%  "
            << G4std::setw(9) << rmsELongit[i] << "%       "
            << setprecision(2)
            << "      0->" << G4std::setw(6) << sup << " radl: "
            << setprecision(4)
            << G4std::setw(7) << MeanELongitCumul[i] << "%  "
            << G4std::setw(7) << rmsELongitCumul[i] << "% "
            << setprecision(2)
            <<G4endl;
   }

  G4cout << G4endl << G4endl << G4endl;

  G4cout << G4endl;
  G4cout << G4std::setw(37) << "SUMMARY" << G4endl;
  G4cout << G4std::setw(42) << "energy deposit : "
         << setprecision(4)
         << G4std::setw(7)  << MeanELongitCumul[nLbin-1] << " % E0 +- "
         << G4std::setw(7)  <<  rmsELongitCumul[nLbin-1] << " % E0" << G4endl;
  G4cout << G4std::setw(42) << "charged traklen: "
         << setprecision(2)
         << G4std::setw(7)  << MeanChargTrLength << " radl +- "
         << G4std::setw(7)  <<  rmsChargTrLength << " radl" << G4endl;
  G4cout << G4std::setw(42) << "neutral traklen: "
         << G4std::setw(7)  << MeanNeutrTrLength << " radl +- "
         << G4std::setw(7)  <<  rmsNeutrTrLength << " radl" << G4endl;


  G4cout.setf(oldform,G4std::ios::floatfield);
  G4cout.precision(oldprec);

  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
