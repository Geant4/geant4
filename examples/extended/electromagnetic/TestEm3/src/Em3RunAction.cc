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
// $Id: Em3RunAction.cc,v 1.18 2002-12-05 17:11:11 gcosmo Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em3RunAction.hh"
#include "Em3RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4Material.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include "g4std/iomanip"

#ifndef G4NOHIST
 #include "AIDA/IAnalysisFactory.h"
 #include "AIDA/ITreeFactory.h"
 #include "AIDA/ITree.h"
 #include "AIDA/IHistogramFactory.h"
 #include "AIDA/IHistogram1D.h"
 #include "AIDA/IAxis.h"
 #include "AIDA/IAnnotation.h"
 #include "AIDA/ITupleFactory.h"
 #include "AIDA/ITuple.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3RunAction::Em3RunAction(Em3DetectorConstruction* det)
:Detector(det)
{
  runMessenger = new Em3RunActionMessenger(this);   

#ifndef G4NOHIST
   // Creating the analysis factory
 IAnalysisFactory* af = AIDA_createAnalysisFactory();
 
 // Creating the tree factory
 ITreeFactory* tf = af->createTreeFactory();
 
 // Creating a tree mapped to an hbook file.
 tree = tf->create("testem3.paw", false, false, "hbook");

 // Creating a histogram factory, whose histograms will be handled by the tree
 hf   = af->createHistogramFactory(*tree);
 
 for (G4int k=0; k<MaxAbsor; k++) histo[k] = 0;
   
 delete tf;
 delete af;
#endif    
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em3RunAction::~Em3RunAction()
{
  delete runMessenger;
  
#ifndef G4NOHIST
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)

  delete hf;
  delete tree;  
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  //
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();
       
  //initialize cumulative quantities
  //
  for (G4int k=0; k<MaxAbsor; k++)     
      sumEAbs[k]=sum2EAbs[k]=sumLAbs[k]=sum2LAbs[k]=0.;

  //drawing
  // 
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
     
  //histograms
  //
  bookHisto();
  
  //example of print dEdx tables
  //
  ////PrintDedxTables();     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3RunAction::bookHisto()
{
#ifndef G4NOHIST
  // book histograms
  const G4String title  = "Edep/Ebeam in absorber ";
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9","10"};
  G4int nbins=100; G4double vmin=0., vmax=1.;
  G4int NbOfAbsor = Detector->GetNbOfAbsor();
  for (G4int k=0; k<NbOfAbsor; ++k)
     {
      if (histo[k]==0)
        { histo[k] = hf->create1D(id[k+1],title+id[k],nbins,vmin,vmax);
          G4cout << "bookHisto: " << k << " " << histo[k] << G4endl;
	}  
     }   
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3RunAction::SetHisto(G4int k,G4int nbins,G4double vmin,G4double vmax)
{
#ifndef G4NOHIST
  // (re)book histograms
  const G4String title  = "Edep/Ebeam in absorber ";
  const G4String id[] = {"0","1","2","3","4","5","6","7","8","9","10"};  
  histo[k] = hf->create1D(id[k+1],title+id[k],nbins,vmin,vmax);
  G4cout << "SetHisto: " << k << " " << histo[k] << G4endl;  
#endif   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


void Em3RunAction::EndOfRunAction(const G4Run* aRun)
{
  if (G4VVisManager::GetConcreteInstance())
     G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewer/update");
     
  //compute and print statistic
  //     
  G4int NbOfEvents = aRun->GetNumberOfEvent();  
  G4double norme = 1./NbOfEvents;
  
  G4double MeanEAbs,rmsEAbs,MeanLAbs,rmsLAbs;
  
#ifdef G4USE_STD_NAMESPACE
  G4std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#else 
  G4long mode = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#endif
  G4int  prec = G4cout.precision(2);
    
  G4cout << "\n-------------------------------------------------------------\n"
         << G4std::setw(51) << "total energy dep" 
	 << G4std::setw(30) << "total tracklen \n \n";
	   
  for (G4int k=0; k<Detector->GetNbOfAbsor(); k++)
    {
     MeanEAbs = norme*sumEAbs[k];
      rmsEAbs = norme*sqrt(abs(NbOfEvents*sum2EAbs[k] - sumEAbs[k]*sumEAbs[k]));
  
     MeanLAbs = norme*sumLAbs[k];
      rmsLAbs = norme*sqrt(abs(NbOfEvents*sum2LAbs[k] - sumLAbs[k]*sumLAbs[k]));
  
     //print
     //    
     G4cout
     << " Absorber" << k 
     << " (" << G4std::setw(12) << Detector->GetAbsorMaterial(k)->GetName() 
     << ") :" 
     << G4std::setw( 7) << G4BestUnit(MeanEAbs,"Energy") << " +- "
     << G4std::setw( 5) << G4BestUnit( rmsEAbs,"Energy")
     << G4std::setw(12) << G4BestUnit(MeanLAbs,"Length") << " +- "
     << G4std::setw( 5) << G4BestUnit( rmsLAbs,"Length")
     << G4endl;
    }
    
  G4cout << "\n-------------------------------------------------------------";
  G4cout << G4endl;  
  G4cout.setf(mode,G4std::ios::floatfield);
  G4cout.precision(prec);
    
  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4EnergyLossTables.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em3RunAction::PrintDedxTables()
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
#ifdef G4USE_STD_NAMESPACE
  G4std::ios::fmtflags mode = G4cout.flags();
  G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#else 
  long mode = G4cout.setf(G4std::ios::fixed,G4std::ios::floatfield);
#endif
  G4int  prec = G4cout.precision(3);
     
  G4cout << "\n kinetic energies \n ";
  for (G4int j=0; j<nbin; ++j) {
    G4cout << G4BestUnit(tk[j],"Energy") << "\t";     
    if ((j+1)%ncolumn == 0) G4cout << "\n ";
  }
  G4cout << G4endl; 

  //print the dE/dx tables
  //
  G4cout.setf(G4std::ios::scientific,G4std::ios::floatfield);
      		    
  G4ParticleDefinition* 
  part = G4ParticleTable::GetParticleTable()->FindParticle("mu+");
  
  for (G4int iab=0;iab < Detector->GetNbOfAbsor(); iab++)
     {
      G4Material* mat = Detector->GetAbsorMaterial(iab);
      G4cout << "\nLIST";
      G4cout << "\nC \nC  dE/dx (MeV/cm) for " << part->GetParticleName() 
             << " in " << mat ->GetName() << "\nC";
      G4cout << "\nKINE   (" << part->GetParticleName() << ")"; 	     
      G4cout << "\nMATE   (" << mat ->GetName() << ")";
      G4cout.precision(2);
      G4cout << "\nERAN  " << tkmin/GeV << " (ekmin)\t"
                           << tkmax/GeV << " (ekmax)\t"
			   << nbin      << " (nekbin)";
      G4double cutgam = (G4Gamma::Gamma()->GetEnergyCuts())[mat->GetIndex()];
      if (cutgam < tkmin) cutgam = tkmin; if (cutgam > tkmax) cutgam = tkmax;
      G4double cutele = (G4Electron::Electron()
                          ->GetEnergyCuts())[mat->GetIndex()];
      if (cutele < tkmin) cutele = tkmin; if (cutele > tkmax) cutele = tkmax;
      G4cout << "\nCUTS  " << cutgam/GeV << " (cutgam)\t" 
                           << cutele/GeV << " (cutele)";
      			   
      G4cout.precision(6);            
      G4cout << "\nG4VAL \n ";
      for (G4int l=0;l<nbin; ++l)
         {
           G4double dedx = G4EnergyLossTables::GetDEDX(part,tk[l],mat);
           G4cout << dedx/(MeV/cm) << "\t";
	   if ((l+1)%ncolumn == 0) G4cout << "\n ";
         }
      G4cout << G4endl; 
     }
     
  G4cout.precision(prec);
  G4cout.setf(mode,G4std::ios::floatfield);     
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
