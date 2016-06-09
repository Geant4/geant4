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
// $Id: RunAction.cc,v 1.24 2008/08/22 18:30:27 vnivanch Exp $
// GEANT4 tag $Name: geant4-09-02 $
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
#include "Histo.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction(DetectorConstruction* det, PhysicsList* phys,
                     PrimaryGeneratorAction* kin)
:detector(det), physics(phys), kinematic(kin)
{ 
  tallyEdep = new G4double[MaxTally];
  binLength = offsetX = 0.;
  histo = new Histo();
  histo->setFileName("testem7");
  histo->add1D("1","Edep (MeV/mm) along absorber (mm)", 100, 0, 100);
  histo->add1D("2","Edep (MeV/mm) along absorber zoomed (mm)", 100, 0, 100);
  histo->add1D("3","Projectile range (mm)", 100, 0, 100);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
  delete [] tallyEdep;
  delete histo;  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::FillHisto(G4int ih, G4double x, G4double weight)
{
  histo->fill(ih, x, weight);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* aRun)
{  
  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;
  
  // save Rndm status
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  CLHEP::HepRandom::showEngineStatus();
     
  //initialize projected range, tallies, Ebeam, and book histograms
  //
  nPrimarySteps = 0;
  nRange = 0;
  projRange = projRange2 = 0.;
  edeptot = eniel = 0.;
  for (G4int j=0; j<MaxTally; j++) tallyEdep[j] = 0.;
  kinematic->ResetEbeamCumul();

  // define "1" histogram binning
  length  = detector->GetAbsorSizeX();
  G4double stepMax = physics->GetStepMaxProcess()->GetMaxStep();
  const G4int nbmin = 100;
  G4int nbBins = (G4int)(0.5 + length/stepMax);
  if (nbBins < nbmin) nbBins = nbmin;
  binLength = length/nbBins;
  offsetX   = 0.5*length;
   
  // histogram "1" is defined by the length of the target
  // zoomed histograms are defined by UI command
  histo->setHisto1D(0, nbBins, 0, length, mm);
 
  histo->book();
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
  if(nRange > 0) {
    projRange /= nRange; 
    projRange2 /= nRange;
  }
  G4double rms = projRange2 - projRange*projRange;        
  if (rms>0.) rms = std::sqrt(rms); else rms = 0.;

  G4double nstep = G4double(nPrimarySteps)/G4double(NbofEvents);

  G4cout.precision(6);       
  G4cout << "\n Projected Range= "<< G4BestUnit(projRange,"Length")
         << "   rms= "            << G4BestUnit( rms,"Length")
         << G4endl;
  G4cout << " Mean number of primary steps = "<< nstep << G4endl;

  //compute energy deposition and NIEL
  //
  edeptot /= NbofEvents; 
  G4cout << " Total energy deposit= "<< G4BestUnit(edeptot,"Energy")
         << G4endl;
  eniel /= NbofEvents; 
  G4cout << " NIEL energy deposit = "<< G4BestUnit(eniel,"Energy")
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
      G4cout << " tally " << j << ": \t \t"
             << G4BestUnit(Edep,"Energy") << "\t"
	     << ratio << " % \t"
	     << G4BestUnit(Dose,"Dose")   << G4endl;
    }
    G4cout << "\n---------------------------------------------------------\n";
    G4cout << G4endl; 
  }

  // normalize histogram
  G4double fac = (mm/MeV)/(NbofEvents *  binLength);
  for (G4int j=0; j<3; j++) {histo->scale(j, fac);}
 
  // save and clean histo
  histo->save();
 
  // show Rndm status
  CLHEP::HepRandom::showEngineStatus();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
