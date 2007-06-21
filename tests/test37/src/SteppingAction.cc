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
#include "SteppingAction.hh"

#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "G4Track.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* arun)
:G4UserSteppingAction(),detector(det),runaction(arun) 
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

SteppingAction::~SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4Track* track = aStep->GetTrack();
  G4VPhysicalVolume* volume = track->GetVolume();
  G4double edep = aStep->GetTotalEnergyDeposit();

  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 1"))
  {// sum per layer
    runaction->SumEnergy1(volume->GetCopyNo(),edep);
  // sum per medium  
    runaction -> AddEnergy1(edep);
  }

  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 2"))
  { runaction->SumEnergy2(volume->GetCopyNo(),edep);
    runaction -> AddEnergy2(edep);
  }
  
  if ((edep!=0.)&&(volume -> GetName() == "LayerMedium 3"))
  { runaction->SumEnergy3(volume->GetCopyNo(),edep);
    runaction -> AddEnergy3(edep);
  }
  
  
}

