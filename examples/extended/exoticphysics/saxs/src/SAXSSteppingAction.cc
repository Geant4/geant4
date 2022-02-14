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
/// \file SAXSSteppingAction.cc
/// \brief Definition of the SAXSSteppingAction class
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4SDManager.hh"
#include "G4UImanager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4AnalysisManager.hh"

#include "SAXSSteppingAction.hh"
#include "SAXSEventAction.hh"
#include "SAXSDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4OpticalPhoton.hh"

#include <cmath>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSSteppingAction::SAXSSteppingAction(SAXSEventAction* eventAction):
  G4UserSteppingAction(),
  fEventAction(eventAction),
  fPhantom(0),
  fEventNumber(-1),
  fNSe(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SAXSSteppingAction::~SAXSSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SAXSSteppingAction::UserSteppingAction(const G4Step* step)
{        
  //DetectorConstruction instance
  const SAXSDetectorConstruction* detectorConstruction = 
    static_cast<const SAXSDetectorConstruction*>
    (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        
  //get Phantom volume
  fPhantom = detectorConstruction->GetPhantom();

  //get pre and post step points
  G4StepPoint* preStepPoint = step->GetPreStepPoint();
  G4StepPoint* postStepPoint = step->GetPostStepPoint();

  //get the volume of the current step
  G4LogicalVolume* volume = preStepPoint->GetTouchableHandle()
    ->GetVolume()->GetLogicalVolume();
          
  //get track and particle name
  G4Track* track = step->GetTrack();
  G4int ID = track->GetTrackID();
  G4String partName = track->GetDefinition()->GetParticleName();

  //update the primary particle (event) weight
  G4double PrimaryWeight = 1.;
  if (partName == "gamma" && ID==1) {
    PrimaryWeight = track->GetWeight();
    fEventAction->UpdateEventWeight(PrimaryWeight);
  }
                                         
  //if the particle is not a primary photon stepping in the phantom, exit!
  if (volume != fPhantom || partName != "gamma" || ID!=1) return;
                
  //if it is a new event, reset variable values to 0
  G4int eventNumber = G4RunManager::GetRunManager()
    ->GetCurrentEvent()->GetEventID();
  if (eventNumber != fEventNumber) {
    fEventNumber = eventNumber;
    fNSe = 0;
  }

  //identify the process that occurred
  G4String Process; 
  if (postStepPoint->GetProcessDefinedStep() != NULL) {
    Process = postStepPoint->GetProcessDefinedStep()->GetProcessName();
  } else {
    Process = "UserLimit";
  }
  //G4cout << "Process: " << Process << G4endl;
  G4int ProcIndex = -1;
  if (Process == "Transportation" || Process == "CoupledTransportation" 
      || Process == "UserLimit") {
    ProcIndex = 0;
  }
  if (Process == "Rayl" || Process == "biasWrapper(Rayl)") {
    ProcIndex = 1;
    fNSe++;
    fEventAction->AddNRi();
  }
  if (Process == "compt" || Process == "biasWrapper(compt)") {
    ProcIndex = 2;
    fNSe++;
    fEventAction->AddNCi();
  }
  if (Process == "phot" || Process == "biasWrapper(phot)") ProcIndex = 3;
  if (Process == "conv" || Process == "biasWrapper(conv)") ProcIndex = 4;        
  if (Process == "photonNuclear" ||
      Process == "biasWrapper(photonNuclear)") {
    ProcIndex = 5;
  }
  if (Process == "PowderDiffraction" ||
      Process == "biasWrapper(PowderDiffraction)") {
    ProcIndex = 6;
    fNSe++;
    fEventAction->AddNDi();
  }
  //G4cout << "ProcIndex: " << ProcIndex << G4endl;
        
  //calculate the scattering angle
  G4ThreeVector mom1 = preStepPoint->GetMomentumDirection();
  G4ThreeVector mom2 = postStepPoint->GetMomentumDirection();        
  G4double theta = mom1.angle(mom2); 

  //fill the Scattering ntuple 
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();
  if (theta > 1e-6) { //record only if the angle is not 0
    analysisManager->FillNtupleIColumn(1,0,ProcIndex);
    analysisManager->FillNtupleDColumn(1,1,
        preStepPoint->GetKineticEnergy()/CLHEP::keV);
    analysisManager->FillNtupleDColumn(1,2,theta/CLHEP::deg); 
    analysisManager->FillNtupleDColumn(1,3,track->GetWeight());    
    analysisManager->AddNtupleRow(1);
  }        
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

