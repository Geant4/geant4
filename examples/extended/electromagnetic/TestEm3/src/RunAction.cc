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
// $Id: RunAction.cc,v 1.16 2004-05-25 20:24:11 vnivanch Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "EmAcceptance.hh"

#include "Randomize.hh"
#include "G4ios.hh"
//#include <iomanip>
//#include <iostream>

#ifdef USE_AIDA
 #include "AIDA/AIDA.h"
#endif

#ifdef USE_ROOT
 #include "TFile.h"
 #include "TH1F.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim)
:Detector(det), Primary(prim)
{
  runMessenger = new RunActionMessenger(this);

  fileName = "testem3.paw";
  for (G4int k=0; k<MaxAbsor; k++) {hbins[k] = 0; histoUnit[k] = 1.;}

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetHisto(G4int k,
               G4int nbins, G4double valmin, G4double valmax, G4String unit)
{
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9","10"};
  G4String title = "Edep in absorber " + id[k] + " (" + unit + ")";
  G4double valunit = G4UnitDefinition::GetValueOf(unit);

  hid[k]    = id[k+1];
  htitle[k] = title;
  hbins[k]  = nbins;
  hmin[k]   = valmin/valunit;
  hmax[k]   = valmax/valunit;
  histoUnit[k] = valunit;
  
#ifdef USE_AIDA  
  G4cout << "---->SetHisto: " << title << " ; " << nbins << " bins from "
         << hmin[k] << " " + unit << " to " << hmax[k] << " " + unit  << G4endl;  
#endif

#ifdef USE_ROOT  
  G4cout << "---->SetHisto: " << title << " ; " << nbins << " bins from "
         << hmin[k] << " " + unit << " to " << hmax[k] << " " + unit  << G4endl;  
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  // save Rndm status
  //
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();

  //initialize cumulative quantities
  //

  for (G4int k=0; k<MaxAbsor; k++) {
      sumEAbs[k] = sum2EAbs[k]  = sumLAbs[k] = sum2LAbs[k] =
      sumEleav[k]= sum2Eleav[k] = 0.;
  }

  //histograms
  //
#ifdef USE_AIDA
  // Creating the analysis factory
  std::auto_ptr<AIDA::IAnalysisFactory> af(AIDA_createAnalysisFactory());

  // Creating the tree factory
  std::auto_ptr<AIDA::ITreeFactory> tf(af->createTreeFactory());

  // Creating a tree mapped to an hbook file.
  G4bool readOnly  = false;
  G4bool createNew = true;
  tree = tf->create(fileName, "hbook", readOnly, createNew);
  
  // Creating a histogram factory, whose histograms will be handled by the tree
  hf = af->createHistogramFactory(*tree);

  // histograms
  for (G4int k=0; k<MaxAbsor; k++) {
   if (hbins[k] > 0)
    histo[k] = hf->createHistogram1D(hid[k],htitle[k],hbins[k],hmin[k],hmax[k]);
   else histo[k] = 0;
  }
#endif

#ifdef USE_ROOT
 // Create a ROOT file
 tree = new TFile(fileName,"recreate");
 
 // histograms
 for (G4int k=0; k<MaxAbsor; k++) {
  if (hbins[k] > 0)
   histo[k] = new TH1F(hid[k],htitle[k],hbins[k],hmin[k],hmax[k]);
  else histo[k] = 0;
 }
#endif 

  //example of print dEdx tables
  //
  ////PrintDedxTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs,
                                         G4double Eleav)
{
  //accumulate statistic
  //
  sumEAbs[kAbs]  += EAbs;  sum2EAbs[kAbs]  += EAbs*EAbs;
  sumLAbs[kAbs]  += LAbs;  sum2LAbs[kAbs]  += LAbs*LAbs;
  sumEleav[kAbs] += Eleav; sum2Eleav[kAbs] += Eleav*Eleav;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double norm = 1.0/(G4double)NbOfEvents;
  G4double qnorm = sqrt(norm);
  G4double beamEnergy = Primary->GetParticleGun()->GetParticleEnergy();
  G4double sqbeam = sqrt(beamEnergy/GeV);

  G4double MeanEAbs,MeanEAbs2,rmsEAbs,resolution,rmsres;
  G4double MeanLAbs,MeanLAbs2,rmsLAbs;
  G4double MeanEleav,MeanEleav2,rmsEleav;

  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
//  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4cout << G4endl;
  G4cout << "-----------------------------------------------------------------------------------------------------------------------"<<G4endl;
  G4cout << "           Material          Total Edep               sqrt(E0(GeV))*rmsE/Emean  "
	 << " total tracklen     Eleak from second Abs" << G4endl;
  G4cout << G4endl;

  for (G4int k=0; k<Detector->GetNbOfAbsor(); k++)
    {
     MeanEAbs  = sumEAbs[k]*norm;
     MeanEAbs2 = sum2EAbs[k]*norm;
      rmsEAbs  = sqrt(abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
     resolution= 100.*sqbeam*rmsEAbs/MeanEAbs;
     rmsres    = resolution*qnorm;

     MeanLAbs  = sumLAbs[k]*norm;
     MeanLAbs2 = sum2LAbs[k]/norm;
      rmsLAbs  = sqrt(abs(MeanLAbs2 - MeanLAbs*MeanLAbs));

     MeanEleav = sumEleav[k]*norm;
     MeanEleav2= sum2Eleav[k]*norm;
      rmsEleav = sqrt(abs(MeanEleav2-MeanEleav*MeanEleav));

     //print
     //
     G4String mat = Detector->GetAbsorMaterial(k)->GetName();
     G4cout
       << "Absorber" << k << " ("
       << std::setw(12) << mat << "): "
       << std::setprecision(5)
       << std::setw(8) << G4BestUnit(MeanEAbs,"Energy")
       << " +- "
       << std::setprecision(4)
       << G4BestUnit( rmsEAbs,"Energy")  << "\t"
       << resolution  << " +- " << rmsres<< " %\t"
       << std::setprecision(3)
       << G4BestUnit(MeanLAbs,"Length")  << " +- "
       << G4BestUnit( rmsLAbs,"Length")  << "\t"
       << G4BestUnit(MeanEleav,"Energy") << " +- "
       << G4BestUnit( rmsEleav,"Energy") << "\t"
       << G4endl;
    }

  G4cout << "----------------------------------------------------------------------------------------------------------------------"<<G4endl;
  G4cout << G4endl;
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // Acceptance
  EmAcceptance acc;
  G4bool isStarted = false;
  for (G4int j=0; j<Detector->GetNbOfAbsor(); j++) {

    G4double ltrue = Detector->GetLimitEdep(j);
    if (ltrue < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Sampling Calorimeter",NbOfEvents);
	isStarted = true;
      }
      G4double etrue = Detector->GetAverageEdep(j);
      G4double rtrue = Detector->GetRMSEdep(j);
      MeanEAbs  = sumEAbs[j]*norm;
      MeanEAbs2 = sum2EAbs[j]*norm;
      rmsEAbs   = sqrt(abs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      G4String mat = Detector->GetAbsorMaterial(j)->GetName();
      acc.EmAcceptanceGauss("Edep"+mat,NbOfEvents,MeanEAbs,etrue,rtrue,ltrue);
      acc.EmAcceptanceGauss("Erms"+mat,NbOfEvents,rmsEAbs,rtrue,rtrue,2.0*ltrue);
    }
  }
  if(isStarted) acc.EndOfAcceptance();

  //save histograms and delete factory
  //
#ifdef USE_AIDA
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
  G4cout << "---> Histograms are saved" << G4endl;
  delete hf;
  delete tree;
#endif

#ifdef USE_ROOT
  tree->Write();        // Writing the histograms to the file
  tree->Close();        // and closing the file
  G4cout << "---> Histograms are saved" << G4endl;
  delete tree;
#endif

  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4LossTableManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::PrintDedxTables()
{
  //Print dE/dx tables with binning identical to the Geant3 JMATE bank.
  //The printout is readable as Geant3 ffread data cards (by the program g4mat).
  //
  const G4double tkmin=10*keV, tkmax=10*TeV;
  const G4int nbin=90;
  G4double tk[nbin];

  const G4int ncolumn = 5;

  //compute the kinetic energies
  //
  const G4double dp = log10(tkmax/tkmin)/nbin;
  const G4double dt = pow(10.,dp);
  tk[0] = tkmin;
  for (G4int i=1; i<nbin; ++i) tk[i] = tk[i-1]*dt;

  //print the kinetic energies
  //
  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(3);

  G4cout << "\n kinetic energies \n ";
  for (G4int j=0; j<nbin; ++j) {
    G4cout << G4BestUnit(tk[j],"Energy") << "\t";
    if ((j+1)%ncolumn == 0) G4cout << "\n ";
  }
  G4cout << G4endl;

  //print the dE/dx tables
  //
  G4cout.setf(std::ios::scientific,std::ios::floatfield);

  G4ParticleDefinition*
  part = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
  const G4ProductionCutsTable* theCoupleTable=
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;

  for (G4int iab=0;iab < Detector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = Detector->GetAbsorMaterial(iab);
      G4int index = 0;
      for ( size_t i=0; i<numOfCouples; i++)
        {
          couple = theCoupleTable->GetMaterialCutsCouple(i);
          if(couple->GetMaterial() == mat) {
	    index = i;
	    break;
	  }
	}
      G4cout << "\nLIST";
      G4cout << "\nC \nC  dE/dx (MeV/cm) for " << part->GetParticleName()
             << " in " << mat ->GetName() << "\nC";
      G4cout << "\nKINE   (" << part->GetParticleName() << ")";
      G4cout << "\nMATE   (" << mat ->GetName() << ")";
      G4cout.precision(2);
      G4cout << "\nERAN  " << tkmin/GeV << " (ekmin)\t"
                           << tkmax/GeV << " (ekmax)\t"
			   << nbin      << " (nekbin)";
      G4double cutgam =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4GammaCut)))[index];
      if (cutgam < tkmin) cutgam = tkmin;
      if (cutgam > tkmax) cutgam = tkmax;
      G4double cutele =
         (*(theCoupleTable->GetEnergyCutsVector(idxG4ElectronCut)))[index];
      if (cutele < tkmin) cutele = tkmin;
      if (cutele > tkmax) cutele = tkmax;
      G4cout << "\nCUTS  " << cutgam/GeV << " (cutgam)\t"
                           << cutele/GeV << " (cutele)";

      G4cout.precision(6);
      G4cout << "\nG4VAL \n ";
      for (G4int l=0;l<nbin; ++l)
         {
           G4double dedx = G4LossTableManager::Instance()
	                                       ->GetDEDX(part,tk[l],couple);
           G4cout << dedx/(MeV/cm) << "\t";
	   if ((l+1)%ncolumn == 0) G4cout << "\n ";
         }
      G4cout << G4endl;
     }

  G4cout.precision(prec);
  G4cout.setf(mode,std::ios::floatfield);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
