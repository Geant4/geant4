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
// $Id: RunAction.cc,v 1.2 2003/11/03 16:42:50 maire Exp $
// GEANT4 tag $Name: geant4-06-00 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>

#ifdef G4ANALYSIS_USE
 #include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det)
:Detector(det)
{
  runMessenger = new RunActionMessenger(this);   

#ifdef G4ANALYSIS_USE
   // Creating the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 AIDA::ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create("testem3.paw", "hbook", readOnly, createNew);

 // Creating a histogram factory, whose histograms will be handled by the tree
 hf   = af->createHistogramFactory(*tree);
 
 for (G4int k=0; k<MaxAbsor; k++) {histo[k] = 0; histoUnit[k] = 1.;}
   
 delete tf;
 delete af;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;

#ifdef G4ANALYSIS_USE
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)

  delete hf;
  delete tree;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetHisto(G4int k,
               G4int nbins, G4double valmin, G4double valmax, G4String unit)
{
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9","10"};
  G4String title = "Edep in absorber " + id[k] + " (" + unit + ")";
  G4double valunit = G4UnitDefinition::GetValueOf(unit);
  valmin /= valunit; valmax /= valunit;
  
  G4cout << "---->SetHisto: " << title << " ; " << nbins << " bins from "
         << valmin << " " << unit << " to " << valmax << unit  << G4endl;
	   
#ifdef G4ANALYSIS_USE  
  histoUnit[k] = valunit;
  histo[k] = hf->createHistogram1D(id[k+1],title,nbins,valmin,valmax);
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
  for (G4int k=0; k<MaxAbsor; k++) 
      sumEAbs[k] = sum2EAbs[k]  = sumLAbs[k] = sum2LAbs[k] = 
      sumEleav[k]= sum2Eleav[k] = 0.;

  //drawing
  //
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");


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
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");

  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double norme = 1./NbOfEvents;

  G4double MeanEAbs,rmsEAbs,MeanLAbs,rmsLAbs,MeanEleav,rmsEleav;

  std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(std::ios::fixed,std::ios::floatfield);
  G4int  prec = G4cout.precision(2);

  G4cout << "\n-------------------------------------------------------------\n"
         << std::setw(51) << "total energy dep"
	 << std::setw(30) << "total tracklen"
	 << std::setw(40) << "mean energy leaving per layer \n \n";	 

  for (G4int k=0; k<Detector->GetNbOfAbsor(); k++)
    {
     MeanEAbs = norme*sumEAbs[k];
      rmsEAbs = norme*sqrt(abs(NbOfEvents*sum2EAbs[k] - sumEAbs[k]*sumEAbs[k]));

     MeanLAbs = norme*sumLAbs[k];
      rmsLAbs = norme*sqrt(abs(NbOfEvents*sum2LAbs[k] - sumLAbs[k]*sumLAbs[k]));
      
     MeanEleav = norme*sumEleav[k];
      rmsEleav=norme*sqrt(abs(NbOfEvents*sum2Eleav[k]-sumEleav[k]*sumEleav[k]));

     //print
     //
     G4cout
     << " Absorber" << k
     << " (" << std::setw(12) << Detector->GetAbsorMaterial(k)->GetName()
     << ") :"
     << std::setw( 7) << G4BestUnit(MeanEAbs,"Energy") << " +- "
     << std::setw( 5) << G4BestUnit( rmsEAbs,"Energy")
     << std::setw(12) << G4BestUnit(MeanLAbs,"Length") << " +- "
     << std::setw( 5) << G4BestUnit( rmsLAbs,"Length")
     << std::setw(22) << G4BestUnit(MeanEleav,"Energy") << " +- "
     << std::setw( 5) << G4BestUnit( rmsEleav,"Energy")     
     << G4endl;
    }

  G4cout << "\n-------------------------------------------------------------";
  G4cout << G4endl;
  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4ProductionCutsTable.hh"
#include "G4EnergyLossTables.hh"

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


  for (G4int iab=0;iab < Detector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = Detector->GetAbsorMaterial(iab);
      G4int index = mat->GetIndex();
      const G4MaterialCutsCouple* couple =
            theCoupleTable->GetMaterialCutsCouple(index);
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
           G4double dedx = G4EnergyLossTables::GetDEDX(part,tk[l],couple);
           G4cout << dedx/(MeV/cm) << "\t";
	   if ((l+1)%ncolumn == 0) G4cout << "\n ";
         }
      G4cout << G4endl;
     }

  G4cout.precision(prec);
  G4cout.setf(mode,std::ios::floatfield);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
