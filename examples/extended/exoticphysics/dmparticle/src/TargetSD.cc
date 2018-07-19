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
/// \file exoticphysics/dmparticle/src/TargetSD.cc
/// \brief Implementation of the TargetSD. class
//
// $Id$
//
/////////////////////////////////////////////////
//
// ClassName:  TargetSD
//
// Authors:  01.06.17 V.Ivanchenko 
//
//
///////////////////////////////////////////

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "TargetSD.hh"
#include "Run.hh"
#include "globals.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4Step.hh"
#include "G4String.hh"
#include "G4RunManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TargetSD::TargetSD(const G4String& name)
  : G4VSensitiveDetector(name), fRun(nullptr)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TargetSD::~TargetSD()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TargetSD::Initialize(G4HCofThisEvent*)
{
  fRun = (Run*)G4RunManager::GetRunManager()->GetNonConstCurrentRun();  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TargetSD::ProcessHits( G4Step* aStep, G4TouchableHistory* )
{
  G4double edep = aStep->GetTotalEnergyDeposit();
  /*
  G4String name = aStep->GetTrack()->GetDefinition()->GetParticleName();

  if(name == "ldmhi" || name == "ldmhibar" || name == "ldmphoton")
  {
  G4cout << name << " Edep(MeV)= " << edep/MeV << G4endl; 
  }
  */
  fRun->AddEnergy(edep, aStep); 

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TargetSD::EndOfEvent(G4HCofThisEvent*)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TargetSD::clear()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TargetSD::PrintAll()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


