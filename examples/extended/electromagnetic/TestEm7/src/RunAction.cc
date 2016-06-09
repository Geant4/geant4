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
// $Id: RunAction.cc,v 1.6 2004/03/31 17:09:46 maire Exp $
// GEANT4 tag $Name: geant4-06-02 $
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

#ifdef USE_AIDA
 #include "AIDA/AIDA.h"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
:detector(det), physics(phys), kinematic(kin)
{ 
  tallyEdep = new G4double[MaxTally];
  binLength = offsetX = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
 cleanHisto();
 delete tallyEdep;
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
 
#ifdef USE_AIDA
 // Create the analysis factory
 AIDA::IAnalysisFactory* af = AIDA_createAnalysisFactory();

 // Create the tree factory
 AIDA::ITreeFactory* tf = af->createTreeFactory();

 // Create a tree mapped to an hbook file.
 G4bool readOnly  = false;
 G4bool createNew = true;
 tree = tf->create("testem7.paw", "hbook", readOnly, createNew);

 // Create a histogram factory, whose histograms will be handled by the tree
 AIDA::IHistogramFactory* hf = af->createHistogramFactory(*tree);

 // Create histograms
 histo[0] = hf->createHistogram1D("1","Edep (MeV/mm)",nbBins, 0,length);

 delete hf;
 delete tf;
 delete af;
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
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  HepRandom::showEngineStatus();
     
  //book histograms and initialize tallies
  //
  if (aRun->GetRunID() == 0) {
    bookHisto();
    for (G4int j=0; j<MaxTally; j++) tallyEdep[j] = 0.;
  }  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run*)
{
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

 // show Rndm status
 HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
