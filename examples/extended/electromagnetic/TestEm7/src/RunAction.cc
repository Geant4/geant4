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
// $Id: RunAction.cc,v 1.14 2005/06/01 13:12:13 maire Exp $
// GEANT4 tag $Name: geant4-07-01 $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "RunAction.hh"
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "StepMax.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4ios.hh"

#include "Randomize.hh"

#ifdef G4ANALYSIS_USE
#include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
:detector(det), physics(phys), kinematic(kin), af(0), tree(0)
{ 
  tallyEdep = new G4double[MaxTally];
  binLength = offsetX = 0.;
  histo[0] = 0;
  
#ifdef G4ANALYSIS_USE
 // Creating the analysis factory
 af = AIDA_createAnalysisFactory();
 if(!af) {
   G4cout << "RunAction::RunAction() :" 
          << " problem creating the AIDA analysis factory."
          << G4endl;
 } 	   
#endif  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete tallyEdep;
  
#ifdef G4ANALYSIS_USE
  delete af;  
#endif      
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::bookHisto()
{
  G4double length  = detector->GetAbsorSizeX();
  G4double stepMax = physics->GetStepMaxProcess()->GetMaxStep();
  const G4int nbmin = 100;
  G4int nbBins = (int)(0.5 + length/stepMax);
  if (nbBins < nbmin) nbBins = nbmin;
  binLength = length/nbBins;
  offsetX   = 0.5*length;
 
#ifdef G4ANALYSIS_USE
  if (!af) return;

  // Create a tree mapped to an hbook file.
  G4bool readOnly  = false;
  G4bool createNew = true;
  G4String options = "--noErrors uncompress";
  AIDA::ITreeFactory* tf  = af->createTreeFactory();  
  tree = tf->create("testem7.hbook","hbook", readOnly, createNew, options);
  //tree = tf->create("testem7.root", "root",readOnly, createNew, options);
  //tree = tf->create("testem7.XML" , "XML" ,readOnly, createNew, options);
  delete tf;
  if (!tree) {
    G4cout << "RunAction::bookHisto()" << G4endl;
    return;
  }

  // Create a histogram factory, whose histograms will be handled by the tree
  AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);
  
  // Create histogram
  histo[0] = hf->createHistogram1D("1","Edep (MeV/mm) along absorber (mm)",
             nbBins, 0, length/mm);
	     
  delete hf;	     
  G4cout << "\n----> Histogram Tree opened" << G4endl;
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::cleanHisto()
{
#ifdef G4ANALYSIS_USE
  tree->commit();       // Writing the histograms to the file
  tree->close();        // and closing the tree (and the file)
  delete tree;
  tree = 0;
  
  G4cout << "\n----> Histogram Tree saved" << G4endl;  
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillHisto(G4int ih, G4double x, G4double weight)
{
#ifdef G4ANALYSIS_USE
  if(histo[ih]) histo[ih]->fill(x, weight);
#endif
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();
     
  //initialize projected range, tallies, Ebeam, and book histograms
  //
  projRange = projRange2 = 0.;
  for (G4int j=0; j<MaxTally; j++) tallyEdep[j] = 0.;
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
   
  G4String particle = kinematic->GetParticleGun()->GetParticleDefinition()
                      ->GetParticleName();    
  G4double energy = kinematic->GetParticleGun()->GetParticleEnergy();
  G4cout << "\n The run consists of " << NbofEvents << " "<< particle << " of "
         << G4BestUnit(energy,"Energy") << " through " 
	 << G4BestUnit(detector->GetAbsorSizeX(),"Length") << " of "
	 << material->GetName() << " (density: " 
	 << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;
	 
  //compute projected range and straggling
  //
  projRange /= NbofEvents; projRange2 /= NbofEvents;
  G4double rms = projRange2 - projRange*projRange;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4cout.precision(5);       
  G4cout << "\n projected Range= "<< G4BestUnit(projRange,"Length")
         << "   rms= "            << G4BestUnit( rms,"Length")
         << G4endl;
     
  //print dose in tallies
  //
  G4int tallyNumber = detector->GetTallyNumber();
  if (tallyNumber > 0) {
    G4double tallyMass = detector->GetTallyMass();
    G4double Ebeam = kinematic->GetEbeamCumul();
    G4cout << "\n---------------------------------------------------------\n";
    G4cout << " Cumulated Doses : \tEdep      \tEdep/Ebeam \tDose" << G4endl;
    for (G4int j=0; j<tallyNumber; j++) {
      G4double Edep = tallyEdep[j], ratio = 100*Edep/Ebeam;
      G4double Dose = Edep/tallyMass;
      G4cout << "tally " << j << ": \t \t"
             << G4BestUnit(Edep,"Energy") << "\t"
	     << ratio << " % \t"
	     << G4BestUnit(Dose,"Dose")   << G4endl;
    }
    G4cout << "\n---------------------------------------------------------\n"; 
  }

#ifdef G4ANALYSIS_USE
  // normalize histogram
  G4double fac = (mm/MeV)/(NbofEvents *  binLength);
  histo[0]->scale(fac);
#endif
  
 
  // save and clean histo
  cleanHisto();
 
  // show Rndm status
  HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
