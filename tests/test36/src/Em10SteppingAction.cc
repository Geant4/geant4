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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

#include <iomanip>

#include "Em10DetectorConstruction.hh"
#include "G4PhysicalConstants.hh"
#include "G4VPhysicalVolume.hh"
#include "G4SteppingManager.hh"
#include "G4TrackVector.hh"
#include "Em10SteppingAction.hh"
#include "Em10PrimaryGeneratorAction.hh"
#include "Em10EventAction.hh"
#include "Em10RunAction.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "Em10SteppingMessenger.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10SteppingAction::Em10SteppingAction(Em10DetectorConstruction* DET,
                                     Em10EventAction* EA,
                                     Em10RunAction* RA)
:detector (DET),eventaction (EA),runaction (RA),steppingMessenger(0),
 IDold(-1) ,evnoold(-1)
{
  steppingMessenger = new Em10SteppingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

Em10SteppingAction::~Em10SteppingAction()
{
  delete steppingMessenger ;
 }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void Em10SteppingAction::UserSteppingAction(const G4Step* aStep)
{ 

  G4int evno = eventaction->GetEventno() ; 

  G4Track* track = aStep->GetTrack();
  G4int trackID = track->GetTrackID();
  G4int parentID = track->GetParentID();
  G4bool isGamma = false;
  G4bool isElectron = false;
  G4bool isPositron = false;
  G4String name = track->GetDefinition()->GetParticleName();
  if("gamma" == name) { isGamma = true; }
  else if("e-" == name) { isElectron = true; }
  else if("e+" == name) { isPositron = true; }

  G4VPhysicalVolume* vol1 = 
    aStep->GetPreStepPoint()->GetPhysicalVolume();
  G4VPhysicalVolume* vol2 = 
    aStep->GetPostStepPoint()->GetPhysicalVolume();
  G4bool absorberStep = false;
  if(vol1 && vol1->GetName() == "absorber") { absorberStep = true; }

  G4double TpreStep = aStep->GetPreStepPoint()->GetKineticEnergy();

  IDnow = evno+10000*trackID+100000000*parentID; 
  if(IDnow != IDold) {
    IDold=IDnow;
    if((isElectron || isPositron) && 0 != parentID) {
      runaction->Fillvertexz(aStep->GetTrack()->GetVertexPosition().x());
    }

    if(absorberStep && 0 != parentID) {
      if(isElectron) {
	eventaction->AddCharged();
	eventaction->AddE();
	runaction->FillTsec(TpreStep);
      } else if(isPositron) {
	eventaction->AddCharged();
	eventaction->AddP();
	runaction->FillTsec(TpreStep);
      } else if(isGamma) {
	eventaction->AddNeutral() ;
      }
    }
  }

  if(absorberStep) {
    if(isElectron || isPositron) {
      eventaction->CountStepsCharged();
    } else if(isGamma) {
      eventaction->CountStepsNeutral();
    }
    
    if(0 == parentID && vol1 != vol2) {
      if(track->GetMomentumDirection().x()>0.) {
	eventaction->SetTr();
	G4double Theta = std::acos(track->GetMomentumDirection().x()) ;
	runaction->FillTh(Theta) ;
	runaction->FillTt(track->GetKineticEnergy()) ;
	G4double yend= track->GetPosition().y() ;
	G4double zend= track->GetPosition().z() ;
	G4double rend = std::sqrt(yend*yend+zend*zend) ;
	runaction->FillR(rend);
      } else {       
	eventaction->SetRef();
	G4double Thetaback = std::acos(track->GetMomentumDirection().x()) ;
	Thetaback -= 0.5*pi ;
	runaction->FillThBack(Thetaback) ;
	runaction->FillTb(track->GetKineticEnergy()) ;
      }
    }
    if (isGamma &&  vol1 != vol2) {
      runaction->FillGammaSpectrum(track->GetKineticEnergy());
    }
  }   
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

