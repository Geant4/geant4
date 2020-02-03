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
/// \file field/field05/src/F05SteppingAction.cc
/// \brief Implementation of the F05SteppingAction class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "F05SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4ProcessTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::F05SteppingAction(void)
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

F05SteppingAction::~F05SteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void F05SteppingAction::UserSteppingAction(const G4Step* aStep)
{

  G4String processName = aStep->GetPostStepPoint()->
                               GetProcessDefinedStep()->GetProcessName();

  if(processName != "DecayWithSpin" ){

     G4Track* aTrack= aStep->GetTrack();

     G4String particleName = aStep->GetTrack()->
                                    GetDefinition()->GetParticleName();

     G4ThreeVector polDir  = aTrack->GetPolarization();
     G4ThreeVector momDir  = aTrack->GetMomentumDirection();

     if (particleName == "mu+") {
        if (momDir * polDir < (1.-1.E-7)) {

           G4double cos_theta = momDir * polDir;
           G4double gTime     = aTrack->GetGlobalTime();

           G4cout << " *** ERROR - WARNING *** " << G4endl;
           G4cout << "processName: " << processName << G4endl;
           G4cout << "particleName " << particleName << G4endl;
           G4cout << "Global Time: " << gTime/ns << "nsec" << G4endl;
           G4cout << "Angle between spin and momentum:" << cos_theta << G4endl;
           G4Exception("SteppingAction::UserSteppingAction","Error",
                       FatalException,
                       "Angle between spin and momentum too large");
        }
     }
  }
}
