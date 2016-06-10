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
/// \file electromagnetic/TestEm10/src/Em10SteppingAction.cc
/// \brief Implementation of the Em10SteppingAction class
//
//
// $Id: Em10SteppingAction.cc 73033 2013-08-15 09:24:45Z gcosmo $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Em10DetectorConstruction.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4Gamma.hh"
#include "Em10SteppingAction.hh"
#include "Em10EventAction.hh"
#include "Em10RunAction.hh"
#include "G4Event.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalConstants.hh"
#include "G4ios.hh"
#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10SteppingAction::Em10SteppingAction(Em10EventAction* EA,
                                       Em10RunAction* RA)
  :G4UserSteppingAction(),eventaction (EA),runaction (RA),
   IDold(-1)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Em10SteppingAction::~Em10SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Em10SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4double Theta,Thetaback,Ttrans,Tback,Tsec,Egamma,yend,zend,rend;

  G4int evno = eventaction->GetEventno();

  const G4Track* track = aStep->GetTrack();
  const G4StepPoint* prePoint = aStep->GetPreStepPoint();
  G4int trackID = track->GetTrackID();
  G4int parentID = track->GetParentID();

  const G4DynamicParticle* dynParticle = track->GetDynamicParticle();
  const G4ParticleDefinition* particle = dynParticle->GetDefinition();
  G4VPhysicalVolume* preVol = prePoint->GetPhysicalVolume();
  G4VPhysicalVolume* postVol = aStep->GetPostStepPoint()->GetPhysicalVolume();

  IDnow = evno+10000*trackID+100000000*parentID;
  if(IDnow != IDold) {
    IDold=IDnow;
    if(trackID > 1 && (particle == G4Electron::Electron() ||
                       particle == G4Positron::Positron() ||
                       particle == G4Gamma::Gamma())) {
      runaction->Fillvertexz(track->GetVertexPosition().z());

      if(preVol->GetName()=="Absorber") {
        if(particle == G4Gamma::Gamma()) {
          eventaction->AddNeutral();
        } else {
          eventaction->AddCharged();
          Tsec = track->GetKineticEnergy();
          Tsec += aStep->GetTotalEnergyDeposit();
          runaction->FillTsec(Tsec);
          if(particle == G4Electron::Electron()) {
            eventaction->AddE();
          } else {
            eventaction->AddP();
          }
        }
      }
    }
  }

  if(preVol->GetName()=="Absorber") {
    if(particle == G4Electron::Electron() ||
       particle == G4Positron::Positron()) {
      eventaction->CountStepsCharged();

    } else if(particle == G4Gamma::Gamma()) {
      eventaction->CountStepsNeutral();
    }

    if(prePoint->GetStepStatus() == fGeomBoundary &&
       preVol != postVol) {

      if(trackID == 1) {
        if(track->GetMomentumDirection().z()>0.) {

          eventaction->SetTr();
          Theta = std::acos(track->GetMomentumDirection().z());
          runaction->FillTh(Theta);
          Ttrans = track->GetKineticEnergy();
          runaction->FillTt(Ttrans);
          yend= aStep->GetTrack()->GetPosition().y();
          zend= aStep->GetTrack()->GetPosition().x();
          rend = std::sqrt(yend*yend+zend*zend);
          runaction->FillR(rend);

        } else {
          eventaction->SetRef();
          Thetaback = std::acos(aStep->GetTrack()->GetMomentumDirection().z());
          Thetaback -= 0.5*pi;
          runaction->FillThBack(Thetaback);
          Tback  = aStep->GetTrack()->GetKineticEnergy();
          runaction->FillTb(Tback);
        }
      }
      if(track->GetMomentumDirection().z()>0. &&
         particle == G4Gamma::Gamma()) {
 
        Egamma = aStep->GetTrack()->GetKineticEnergy();
        runaction->FillGammaSpectrum(Egamma);
      }
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
