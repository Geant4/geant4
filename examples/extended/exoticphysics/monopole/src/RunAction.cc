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
// $Id: RunAction.cc,v 1.5 2010-11-09 21:33:25 asaim Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "QGSP.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "Randomize.hh"
#include "G4EmCalculator.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorAction* kin)
  :detector(det), kinematic(kin), af(0), tree(0)
{ 
  verboseLevel = 1;
  binLength = offsetX = 0.;
  histo[0] = 0;
  tree = 0;
  af   = 0;  
#ifdef G4ANALYSIS_USE
  // Creating the analysis factory
  af = AIDA_createAnalysisFactory();
  ftype   = "root";
  fname   = "monopole";
#endif

  // create commands for interactive definition of the detector  
  runActionMessenger = new RunActionMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
#ifdef G4ANALYSIS_USE
  delete af;  
#endif      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
  G4double length  = detector->GetAbsorSizeX();
  if(!binLength) { binLength = 5 * mm; }
  if(binLength > detector->GetMaxStepSize()) { 
    binLength = detector->GetMaxStepSize();
  }
  offsetX   = 0.5 * length;
 
#ifdef G4ANALYSIS_USE
  if(GetVerbose() > 0) { G4cout << "\n----> Histogram Tree opened" << G4endl; }

  G4int nbBins = (G4int)(0.5 + length / binLength);

  // Create the tree factory
  AIDA::ITreeFactory* tf  = af->createTreeFactory();

  // Create a tree mapped to an hbook file.
  G4bool readOnly  = false;
  G4bool createNew = true;
  //G4String ftype   = "hbook";
  //G4String fname   = "monopole";
	G4String fName   = fname;
  fName += ".";
  fName += ftype;
  G4String option  = "";
  tree = tf->create(fName,ftype, readOnly, createNew, option);

  // Create a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

  // Create histograms
  histo[0] = hf->createHistogram1D("1","Edep (MeV/mm) along absorber (mm)", nbBins, 0, length);
  histo[1] = hf->createHistogram1D("2","DEDX (MeV/mm) of proton", 100, -3., 7.);
  histo[2] = hf->createHistogram1D("3","DEDX (MeV/mm) of monopole", 100, -3., 7.);
  histo[3] = hf->createHistogram1D("4","Range(mm) of proton", 100, -3., 7.);
  histo[4] = hf->createHistogram1D("5","Range(mm) of monopole", 100, -3., 7.);

  delete tf;
  delete hf;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::saveHisto()
{
#ifdef G4ANALYSIS_USE
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
  delete tree;
  tree = 0;
  if(GetVerbose() > 0) G4cout << "\n----> Histogram Tree saved" << G4endl;  
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::SetBinSize(G4double size)
{
  binLength =  size;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillHisto(G4int ih, G4double x, G4double weight)
{
  if(GetVerbose() > 1) {
    G4cout << "FillHisto " << ih << "  x=" << x << " weight= " << weight 
	   << G4endl;
  }
#ifdef G4ANALYSIS_USE
  if(histo[ih]) histo[ih]->fill(x, weight);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  if(GetVerbose() > 0) G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
     
  //initialize projected range, tallies, Ebeam, and book histograms
  projRange = projRange2 = 0.;
  kinematic->ResetEbeamCumul();
  bookHisto();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int NbofEvents = aRun->GetNumberOfEvent();
  if (NbofEvents == 0) return;

  //run conditions
  //  
  G4Material* material = detector->GetAbsorMaterial();
  G4double density = material->GetDensity();
  const G4ParticleDefinition* part = 
    kinematic->GetParticleGun()->GetParticleDefinition();
  G4String particle = part->GetParticleName();    
  G4double energy = kinematic->GetParticleGun()->GetParticleEnergy();

  if(GetVerbose() > 0){
    G4cout << "\n The run consists of " << NbofEvents << " "<< particle << " of "
           << G4BestUnit(energy,"Energy") << " through " 
	   << G4BestUnit(detector->GetAbsorSizeX(),"Length") << " of "
	   << material->GetName() << " (density: " 
	   << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
  };
	 
  //compute projected range and straggling

  projRange /= NbofEvents; projRange2 /= NbofEvents;
  G4double rms = projRange2 - projRange*projRange;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  if(GetVerbose() > 0){
    G4cout.precision(5);       
    G4cout << "\n projected Range= " << G4BestUnit(projRange, "Length")
           << "   rms= "             << G4BestUnit(rms, "Length")
           << G4endl;
  };

  G4double ekin[100], dedxproton[100], dedxmp[100];
  G4EmCalculator calc;
  calc.SetVerbose(0);
  G4int i;
  for(i = 0; i < 100; ++i) {
    ekin[i] = std::pow(10., 0.1*G4double(i)) * keV;
    dedxproton[i] = calc.ComputeElectronicDEDX(ekin[i], "proton",  material->GetName());
    dedxmp[i] = calc.ComputeElectronicDEDX(ekin[i], "monopole",  material->GetName());
  }

  if(GetVerbose() > 0){
    G4cout << "### Stopping Powers" << G4endl;
    for(i=0; i<100; i++) {
      G4cout << " E(MeV)= " << ekin[i] << "  dedxp= " << dedxproton[i]
	     << " dedxmp= " << dedxmp[i]
	     << G4endl;
    }
  };
  G4cout << "### End of stopping power table" << G4endl;
#ifdef G4ANALYSIS_USE
  // normalize histogram
  G4double fac = (mm/MeV) / (NbofEvents * binLength);
  histo[0]->scale(fac);

	G4String matName = detector->GetAbsorMaterial()->GetName();
	if(GetVerbose() > 0){
		G4cout << "Range table for " << matName << G4endl;
	};

  for(i=0; i<100; ++i) {
    G4double e = std::log10(ekin[i] / MeV) + 0.05;
    histo[1]->fill(e, dedxproton[i]);
    histo[2]->fill(e, dedxmp[i]);
    histo[3]->fill(e, std::log10(calc.GetRange(ekin[i], "proton", matName) / mm));
    histo[4]->fill(e, std::log10(calc.GetRange(ekin[i], "monopole", matName) / mm));
  }
  
#endif
 
  // save and clean histo
  saveHisto();
 
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
