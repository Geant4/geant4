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
// $Id: RunAction.cc,v 1.23 2004/12/02 16:13:47 vnivanch Exp $
// GEANT4 tag $Name: geant4-07-00-cand-03 $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"

#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"
#include "HistoManager.hh"
#include "EmAcceptance.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* prim,
                     HistoManager* hist)
:Detector(det), Primary(prim), histoManager(hist)
{
  runMessenger = new RunActionMessenger(this);

  for (G4int k=0; k<MaxAbsor; k++) { edeptrue[k] = rmstrue[k] = 1.;
                                    limittrue[k] = DBL_MAX;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete runMessenger;
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
      sumEAbs[k] = sum2EAbs[k]  = sumLAbs[k] = sum2LAbs[k] = 0.;
  }
  
  //histograms
  //
  histoManager->book();
   
  //example of print dEdx tables
  //
  ////PrintDedxTables();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::fillPerEvent(G4int kAbs, G4double EAbs, G4double LAbs)
{
  //accumulate statistic
  //
  sumEAbs[kAbs]  += EAbs;  sum2EAbs[kAbs]  += EAbs*EAbs;
  sumLAbs[kAbs]  += LAbs;  sum2LAbs[kAbs]  += LAbs*LAbs;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void RunAction::EndOfRunAction(const G4Run* aRun)
{
  //compute and print statistic
  //
  G4int NbOfEvents = aRun->GetNumberOfEvent();
  G4double  norm = 1.0/NbOfEvents;
  G4double qnorm = std::sqrt(norm);
  G4double beamEnergy = Primary->GetParticleGun()->GetParticleEnergy();
  G4double sqbeam = std::sqrt(beamEnergy/GeV);

  G4double MeanEAbs,MeanEAbs2,rmsEAbs,resolution,rmsres;
  G4double MeanLAbs,MeanLAbs2,rmsLAbs;

  std::ios::fmtflags mode = G4cout.flags();
  G4int  prec = G4cout.precision(2);
  G4cout << "\n------------------------------------------------------------\n";
  G4cout << std::setw( 8) << "material"
         << std::setw(23) << "Total Edep"
	 << std::setw(27) << "sqrt(E0(GeV))*rmsE/Emean"
	 << std::setw(23) << "total tracklen \n \n";

  for (G4int k=1; k<=Detector->GetNbOfAbsor(); k++)
    {
     MeanEAbs  = sumEAbs[k]*norm;
     MeanEAbs2 = sum2EAbs[k]*norm;
      rmsEAbs  = std::sqrt(std::fabs(MeanEAbs2 - MeanEAbs*MeanEAbs));
     resolution= 100.*sqbeam*rmsEAbs/MeanEAbs;
     rmsres    = resolution*qnorm;

     MeanLAbs  = sumLAbs[k]*norm;
     MeanLAbs2 = sum2LAbs[k]*norm;
      rmsLAbs  = std::sqrt(std::fabs(MeanLAbs2 - MeanLAbs*MeanLAbs));

     //print
     //
     G4cout
       << std::setw(10) << Detector->GetAbsorMaterial(k)->GetName() << ": "
       << std::setprecision(5)
       << std::setw(6) << G4BestUnit(MeanEAbs,"Energy") << " +- "
       << std::setprecision(4)
       << std::setw(5) << G4BestUnit( rmsEAbs,"Energy")  
       << std::setw(10) << resolution  << " +- " 
       << std::setw(5) << rmsres << " %"
       << std::setprecision(3)
       << std::setw(10) << G4BestUnit(MeanLAbs,"Length")  << " +- "
       << std::setw(4) << G4BestUnit( rmsLAbs,"Length") 
       << G4endl;
    }
  G4cout << "\n------------------------------------------------------------\n";

  G4cout.setf(mode,std::ios::floatfield);
  G4cout.precision(prec);

  // Acceptance
  EmAcceptance acc;
  G4bool isStarted = false;
  for (G4int j=1; j<=Detector->GetNbOfAbsor(); j++) {
    if (limittrue[j] < DBL_MAX) {
      if (!isStarted) {
        acc.BeginOfAcceptance("Sampling Calorimeter",NbOfEvents);
	isStarted = true;
      }
      MeanEAbs  = sumEAbs[j]*norm;
      MeanEAbs2 = sum2EAbs[j]*norm;
      rmsEAbs   = std::sqrt(std::fabs(MeanEAbs2 - MeanEAbs*MeanEAbs));
      G4String mat = Detector->GetAbsorMaterial(j)->GetName();
      acc.EmAcceptanceGauss("Edep"+mat, NbOfEvents, MeanEAbs,
                             edeptrue[j], rmstrue[j], limittrue[j]);
      acc.EmAcceptanceGauss("Erms"+mat, NbOfEvents, rmsEAbs,
                             rmstrue[j], rmstrue[j], 2.0*limittrue[j]);
    }
  }
  if(isStarted) acc.EndOfAcceptance();

  //normalize histograms
  //
  G4double fac = 1./NbOfEvents;
  for (G4int ih = MaxAbsor+1; ih < MaxHisto; ih++) {
     histoManager->Normalize(ih,fac/MeV);
  }
  
  //save histograms   
  histoManager->save();

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
  const G4double dp = std::log10(tkmax/tkmin)/nbin;
  const G4double dt = std::pow(10.,dp);
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
  
  G4ProductionCutsTable* theCoupleTable =
        G4ProductionCutsTable::GetProductionCutsTable();
  size_t numOfCouples = theCoupleTable->GetTableSize();
  const G4MaterialCutsCouple* couple = 0;

  for (G4int iab=0;iab < Detector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = Detector->GetAbsorMaterial(iab);
      G4int index = 0;
      for (size_t i=0; i<numOfCouples; i++) {
         couple = theCoupleTable->GetMaterialCutsCouple(i);
         if (couple->GetMaterial() == mat) {index = i; break;}
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

void RunAction::SetEdepAndRMS(G4int i, G4ThreeVector Value)
{
  if (i>=0 && i<MaxAbsor) {
    edeptrue [i] = Value(0);
    rmstrue  [i] = Value(1);
    limittrue[i] = Value(2);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
