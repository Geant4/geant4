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
/// \file runAndEvent/RE04/src/RE04SteppingAction.cc
/// \brief Implementation of the RE04SteppingAction class
//
//
#include "RE04SteppingAction.hh"

#include "G4Step.hh"
#include "G4Track.hh"
#include "G4Material.hh"
#include "G4StepPoint.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04SteppingAction::RE04SteppingAction()
 : G4UserSteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
RE04SteppingAction::~RE04SteppingAction()
{;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void RE04SteppingAction::UserSteppingAction(const G4Step* aStep)
{
  const G4Track* trk = aStep->GetTrack();
  G4cout << "PreStep " << aStep->GetPreStepPoint()->GetPosition()
         << " -- PostStep " << aStep->GetPostStepPoint()->GetPosition()
         << " >>> Material " << trk->GetMaterial()->GetName() << G4endl;
}

